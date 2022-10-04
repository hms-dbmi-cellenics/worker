mock_req <- function() {
  req <- list(body = list(genes = list("MS4A1", "CD79B")))
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
    c("order", "stats", "rawExpression", "truncatedExpression", "zScore")
  )
  expect_equal(length(res$order), 2)
  expect_equal(length(res$stats), 2)
  expect_equal(res$order, req$body$genes)
  expect_equal(names(res$stats), unlist(req$body$genes))
})

test_that("Expression matrices are correctly formatted for mathJS", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression,
    trunc = res$truncatedExpression,
    z = res$zScore
  )
  # expression matrices contain 4 attributes
  expect_equal(unlist(unique(lapply(exp, length))), 4)

  # expression matrices contain correctly named attributes
  expect_equal(
    unlist(unique(lapply(exp, names))),
    c("values", "index", "ptr", "size")
  )

  # no NA values
  expect_true(!all(unlist(lapply(res$rawExpression, is.na))))
  expect_true(!all(unlist(lapply(res$truncatedExpression, is.na))))
  expect_true(!all(unlist(lapply(res$zScore, is.na))))
})


test_that("Expression task returns appropriate number of cells.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression,
    trunc = res$truncatedExpression,
    z = res$zScore
  )

  expected_cells <- max(data$cells_id) + 1

  # assuming all matrices have to be the same size (they have to be)
  result_size <- unique(unlist(lapply(exp, function(x) x$size[[1]])))

  expect_equal(result_size, expected_cells)
})

test_that("summaryStats calculates correct summary stats.", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )


  expression_values <- getExpressionValues(data, gene_subset)
  res <- summaryStats(expression_values)

  expect_equal(names(res), unlist(req$body$genes))
  expect_equal(
    unique(unlist(lapply(res, names))),
    c("rawMean", "rawStdev", "truncatedMin", "truncatedMax")
  )

  expect_equal(
    res$MS4A1$rawMean,
    mean(expression_values$rawExpression$MS4A1, na.rm = TRUE)
  )
  expect_equal(
    res$MS4A1$rawStdev,
    sd(expression_values$rawExpression$MS4A1, na.rm = TRUE)
  )
  expect_equal(
    res$MS4A1$truncatedMin,
    min(expression_values$truncatedExpression$MS4A1)
  )
  expect_equal(
    res$MS4A1$truncatedMax,
    max(expression_values$truncatedExpression$MS4A1)
  )

  expect_equal(
    res$CD79B$rawMean,
    mean(expression_values$rawExpression$CD79B, na.rm = TRUE)
  )
  expect_equal(
    res$CD79B$rawStdev,
    sd(expression_values$rawExpression$CD79B, na.rm = TRUE)
  )
  expect_equal(
    res$CD79B$truncatedMin,
    min(expression_values$truncatedExpression$CD79B)
  )
  expect_equal(
    res$CD79B$truncatedMax,
    max(expression_values$truncatedExpression$CD79B)
  )
})


test_that("Expression task does not return gene that does not exist", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- append(req$body$genes, "aaa")

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression,
    trunc = res$truncatedExpression,
    z = res$zScore
  )

  # non-existent gene is not present in the result
  expect_true(!("aaa" %in% unlist(res$order)))
  expect_true(!("aaa" %in% names(res$stats)))

  expect_equal(length(res$order), 2)
  expect_equal(length(res$stats), 2)

  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[2]]))), 2)
})

test_that("Expression task works with one gene", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("MS4A1")

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression,
    trunc = res$truncatedExpression,
    z = res$zScore
  )
  expect_equal(res$order, req$body$genes)
  expect_equal(length(res$order), 1)
  expect_equal(length(res$stats), 1)

  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[2]]))), 1)
})

test_that("If max truncated expression value is 0, finds a non-zero value", {
  data <- mock_scdata()
  req <- mock_req()

  data@assays$RNA@data["MS4A1", ] <- 0
  data@assays$RNA@data["MS4A1", 1] <- 5

  res <- runExpression(req, data)
  expect_false(all(res$truncatedExpression$values == 0))
})

test_that("runExpression throws an error if request only non existing genes", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("blah")

  expect_error(runExpression(req, data), "Gene\\(s\\): blah not found!")
})


test_that("truncateExpression truncates correctly", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  res <- getExpressionValues(data, gene_subset)

  # not pretty, but we would expect that the max raw expression value per gene
  # to be gte than the adjusted max
  max_raw <- apply(res$rawExpression, 2, max)
  max_adj <- apply(res$truncatedExpression, 2, max)
  expect_true(all(max_raw >= max_adj))

  # check that the raw and truncated matrices are different
  expect_false(isTRUE(all.equal(res$rawExpression, res$truncatedExpression)))

  # truncation effectively makes values equal. So there are less unique values
  # in the truncated table
  expect_gt(
    length(unique(unlist(res$rawExpression))),
    length(unique(unlist(res$truncatedExpression)))
  )
})


test_that("scaleExpression correctly calculates zScore", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  res <- getExpressionValues(data, gene_subset)
  cols <- colnames(res$rawExpression)
  zScore <- data.table::copy(res$rawExpression)

  # calculate zScore in an alternative way
  calculate_zscore <- function(x) {
    as.vector(scale(x))
  }

  zScore[, (cols) := lapply(.SD, calculate_zscore), .SDcols = cols]

  expect_equal(res$zScore, zScore)
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
  original_data <- Matrix::t(data@assays$RNA@data[unlist(req$body$genes), , drop = FALSE])
  original_data <- as.data.table(original_data)

  res <- getExpressionValues(data, gene_subset)

  expect_equal(
    res$rawExpression,
    original_data
  )
})


test_that("order of cells in the completed matrix is correct", {
  data <- mock_scdata()

  cell_ids <- c(0, 5, 10, 30, 79)
  data <- subsetIds(data, cell_ids)

  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  original_data <- Matrix::t(data@assays$RNA@data[unlist(req$body$genes), ,
    drop = FALSE
  ])
  original_data <- data.table::as.data.table(original_data)

  expression_values <- getExpressionValues(data, gene_subset)

  # check that the subsetted data is in the same order as the original data
  # by checking the values directly
  expect_equal(
    colnames(expression_values$rawExpression),
    colnames(original_data)
  )

  expect_equal(expression_values$rawExpression, original_data)

  filled_expression <- completeExpression(
    expression_values$rawExpression,
    data@meta.data$cells_id
  )

  # check that the expression values for each cell ID are located at that cell
  # index (+1 because cell ids are 0-indexed)
  expect_equal(filled_expression[cell_ids + 1, ], original_data)
  # check that all other values are NA
  expect_true(all(is.na(filled_expression[-(cell_ids + 1), ])))
})

test_that("toSparseJson returns values as lists", {
  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  res_long <- getExpressionValues(data, gene_subset)
  res_long_sparse <- lapply(res_long, sparsify)
  res_long_json <- lapply(res_long_sparse, toSparseJson)

  # test that every element in the res_long_json is composed of lists
  lapply(res_long_json, \(x) {
    lapply(x, expect_type, "list")
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

  data_short <- subsetIds(data, 0)
  gene_subset <- gene_subset[1, ]

  res_short <- getExpressionValues(data_short, gene_subset)
  res_short_sparse <- lapply(res_short, sparsify)
  res_short_json <- lapply(res_short_sparse, toSparseJson)

  # test that every element in the res_long_json is composed of lists
  lapply(res_short_json, \(x) {
    lapply(x, expect_type, "list")
  })
})
