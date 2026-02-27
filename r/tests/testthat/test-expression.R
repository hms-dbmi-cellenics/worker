mock_req <- function() {
  req <- list(
    body = list(
      genes = list("MS4A1", "CD79B"),
      downsampled = FALSE
    )
  )
  return(req)
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )
  return(pbmc_small)
}

test_that("Expression task returns appropriate number and names of genes.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  expect_equal(
    names(res),
    c("orderedGeneNames", "stats", "rawExpression")
  )
  expect_equal(length(res$orderedGeneNames), 2)
  expect_equal(length(res$stats$rawMean), 2)
  expect_equal(length(res$stats$rawStdev), 2)
  expect_equal(length(res$stats$truncatedMin), 2)
  expect_equal(length(res$stats$truncatedMax), 2)
  expect_equal(as.list(res$orderedGeneNames), req$body$genes)
})

test_that("Expression task works correctly with downsampled = TRUE.", {
  data <- mock_scdata()
  req <- mock_req()

  req$body$downsampled = TRUE
  req$body$downsampleSettings = list(
    selectedCellSet = "louvain",
    groupedTracks = list(
      "louvain",
      "sample"
    ),
    selectedPoints = "All",
    hiddenCellSets = list()
  )
  req$body$cellIds <- c(1,2,3,4,5)

  res <- runExpression(req, data)

  expect_equal(
    names(res),
    c("orderedGeneNames", "stats", "rawExpression")
  )

  expect_equal(res$orderedGeneNames, c("MS4A1", "CD79B"))
  expect_equal(res$stats$rawMean, c(0.7890259, 1.3545382))
  expect_snapshot(res)
})

test_that("Expression task keeps order regardless of the request received.", {
  data <- mock_scdata()
  req <- mock_req()

  # Change order so that we can make sure this doesnt affect order
  rev_req <- req
  rev_req$body_genes<- rev(req$body$genes)

  res <- runExpression(req, data)
  rev_res <- runExpression(rev_req, data)

  expect_equal(res, rev_res)

  expect_equal(
    names(res),
    c("orderedGeneNames", "stats", "rawExpression")
  )

  expect_equal(res$orderedGeneNames, c("MS4A1", "CD79B"))
  expect_equal(res$stats$rawMean, c(0.7890259, 1.3545382))
  expect_snapshot({
    res
    rev_res})
})

test_that("Expression matrices are correctly formatted for mathJS", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression
  )
  # expression matrix contains 4 attributes
  expect_equal(unlist(unique(lapply(exp, length))), 4)

  # expression matrix contains correctly named attributes
  expect_equal(
    unlist(unique(lapply(exp, names))),
    c("values", "index", "ptr", "size")
  )

  # no NA values
  expect_true(!all(unlist(lapply(res$rawExpression, is.na))))
})


test_that("Expression task returns appropriate number of cells.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression
  )

  expected_cells <- max(data$cells_id) + 1

  # assuming expression matrix has to be the correct size
  result_size <- unique(unlist(lapply(exp, function(x) x$size[[1]])))

  expect_equal(result_size, expected_cells)
})

test_that("getStats works well", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  expression_values <- getExpressionValues(data, gene_subset)
  res <- getStats(expression_values)

  expect_equal(
    names(res),
    c("rawMean", "rawStdev", "truncatedMin", "truncatedMax")
  )

  expect_equal(
    res$rawMean[1],
    mean(expression_values[, 1], na.rm = TRUE)
  )
  expect_equal(
    res$rawStdev[1],
    sd(as.numeric(expression_values[, 1]), na.rm = TRUE)
  )
  expect_equal(
    res$truncatedMin[1],
    min(expression_values[, 1], na.rm = TRUE)
  )
  expect_equal(
    res$truncatedMax[1],
    getQuantileCap(expression_values[, 1, drop = FALSE], 0.95)
  )

  expect_equal(
    res$rawMean[2],
    mean(expression_values[, 2], na.rm = TRUE)
  )
  expect_equal(
    res$rawStdev[2],
    sd(as.numeric(expression_values[, 2]), na.rm = TRUE)
  )
  expect_equal(
    res$truncatedMin[2],
    min(expression_values[, 2], na.rm = TRUE)
  )
  expect_equal(
    res$truncatedMax[2],
    getQuantileCap(expression_values[, 2, drop = FALSE], 0.95)
  )
})


test_that("Expression task does not return gene that does not exist", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- append(req$body$genes, "aaa")

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression
  )

  # non-existent gene is not present in the result
  expect_true(!("aaa" %in% unlist(res$orderedGeneNames)))
  expect_true(!("aaa" %in% names(res$stats)))

  expect_equal(length(res$orderedGeneNames), 2)
  expect_equal(length(res$stats$rawMean), 2)
  expect_equal(length(res$stats$rawStdev), 2)
  expect_equal(length(res$stats$truncatedMin), 2)
  expect_equal(length(res$stats$truncatedMax), 2)

  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[2]]))), 2)
})

test_that("Expression task works with one gene", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("MS4A1")

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression
  )
  expect_equal(res$orderedGeneNames, req$body$genes)
  expect_equal(length(res$orderedGeneNames), 1)
  expect_equal(length(res$stats$rawMean), 1)
  expect_equal(length(res$stats$rawStdev), 1)
  expect_equal(length(res$stats$truncatedMin), 1)
  expect_equal(length(res$stats$truncatedMax), 1)


  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[2]]))), 1)
})

test_that("runExpression throws an error if request only non existing genes", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("blah")

  expect_error(runExpression(req, data), "Gene\\(s\\): blah not found!")
})


test_that("raw expression values are the same as in the original data", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  # manually extract original data, and convert to data.table
  original_data <- Matrix::t(data@assays$RNA$data[unlist(req$body$genes), , drop = FALSE])
  original_data <- as.data.table(original_data)

  res <- getExpressionValues(data, gene_subset)

  # Convert sparse matrix to data.table for comparison
  res_dt <- as.data.table(as.matrix(res))
  
  expect_equal(
    res_dt,
    original_data
  )
})


test_that("toSparseJson returns vectors of lenght > 1 as vectors", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  res_long <- getExpressionValues(data, gene_subset)
  res_long_json <- toSparseJson(res_long)

  # test that every element in the res_long_json is numeric or list
  lapply(res_long_json, \(x){
    expect_true(is.numeric(x) || is.list(x))
  })
})


test_that("toSparseJson returns single value arrays as lists", {

  # RJSONIO::toJSON converts single values to JS scalars. Which breaks when the
  # UI expects single value arrays. Converting all vectors to list using as.list
  # fixes this issue.

  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  # Use 2-3 cells and 1 gene to test formatting
  cell_ids <- c(0, 5)
  data_short <- subsetIds(data, cell_ids)
  gene_subset <- gene_subset[1, ]

  res_short <- getExpressionValues(data_short, gene_subset)
  res_short_json <- toSparseJson(res_short)

  # test that elements are properly formatted
  lapply(res_short_json, \(x) {
    expect_true(is.numeric(x) || is.list(x))
  })
})

