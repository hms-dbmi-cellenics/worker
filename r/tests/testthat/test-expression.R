mock_req <- function() {
  req <- list(
    body = list(
      genes = list("MS4A1", "CD79B"),
      downsampled = FALSE
    )
  )
  return(req)
}

test_that("Expression task works with bpcells.", {
  data <- mock_scdata(use_bpcells = TRUE)
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

test_that("Expression task returns correct values.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

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
  original_data <- Matrix::t(
    data@assays$RNA$data[gene_subset$input, , drop = FALSE]
  )
  colnames(original_data) <- gene_subset$name
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


test_that("getStats produces consistent results with master version", {
  # Define the old getStats from master for comparison
  getStatsOld <- function(data) {
    stats_unsafe <- list(
      rawMean = unname(colMeans(data$rawExpression, na.rm = TRUE)),
      rawStdev = unname(apply(data$rawExpression, 2, sd, na.rm = TRUE)),
      truncatedMin = unname(apply(data$truncatedExpression, 2, min, na.rm = TRUE)),
      truncatedMax = unname(apply(data$truncatedExpression, 2, max, na.rm = TRUE))
    )
    return(stats_unsafe)
  }

  data <- mock_scdata()
  req <- mock_req()
  gene_annotations <- data@misc$gene_annotations

  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  expression_values <- getExpressionValues(data, gene_subset)
  
  # Create old-style data structure with truncatedExpression for comparison
  old_style_data <- list(
    rawExpression = as.data.table(as.matrix(expression_values)),
    truncatedExpression = as.data.table(as.matrix(expression_values))
  )
  
  new_results <- getStats(expression_values)
  old_results <- getStatsOld(old_style_data)
  
  # Compare rawMean and rawStdev (should match)
  expect_equal(new_results$rawMean, old_results$rawMean, tolerance = 1e-10)
  expect_equal(new_results$rawStdev, old_results$rawStdev, tolerance = 1e-10)
  
  # truncatedMin should match (both use min of expression values)
  expect_equal(new_results$truncatedMin, old_results$truncatedMin, tolerance = 1e-10)
})

test_that("expandMatrixToCellIDs creates correct dimensions and zero rows", {
  # Create a small sparse matrix with known values (4 rows, 2 columns)
  i <- c(1, 2, 3, 4)  # 1-based row indices
  j <- c(1, 2, 1, 2)  # 1-based column indices
  x <- c(1.0, 2.0, 3.0, 4.0)  # values
  
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(4, 2))
  
  # Cell IDs: 0, 1, 3, 100 (gaps and large index)
  all_cell_ids <- c(0, 1, 3, 100)
  
  result <- expandMatrixToCellIDs(mat, all_cell_ids)
  
  # Should have 101 rows (0-indexed to 100) and 2 columns
  expect_equal(nrow(result), 101)
  expect_equal(ncol(result), 2)
})

test_that("expandMatrixToCellIDs maps values to correct cell_id rows", {
  # Create a sparse matrix with known values
  # Row 1: [1.0, 0]
  # Row 2: [0, 2.0]
  # Row 3: [3.0, 0]
  # Row 4: [0, 4.0]
  i <- c(1, 2, 3, 4)
  j <- c(1, 2, 1, 2)
  x <- c(1.0, 2.0, 3.0, 4.0)
  
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(4, 2))
  
  # cell_ids: 0, 1, 3, 100
  # Result should be:
  # Row 0 (cell_id 0, from expr row 1): [1.0, 0]
  # Row 1 (cell_id 1, from expr row 2): [0, 2.0]
  # Row 2 (cell_id 2, no expr): [0, 0] (implicit in sparse matrix)
  # Row 3 (cell_id 3, from expr row 3): [3.0, 0]
  # Rows 4-99 (cell_ids 4-99, no expr): [0, 0]
  # Row 100 (cell_id 100, from expr row 4): [0, 4.0]
  
  all_cell_ids <- c(0, 1, 3, 100)
  result <- expandMatrixToCellIDs(mat, all_cell_ids)
  
  # Convert to dense for easier testing of specific rows
  result_dense <- as.matrix(result)
  
  # Check row 0 (cell_id 0, from expression row 1)
  expect_equal(result_dense[1, ], c(1.0, 0))
  
  # Check row 1 (cell_id 1, from expression row 2)
  expect_equal(result_dense[2, ], c(0, 2.0))
  
  # Check row 2 (cell_id 2, should be all zeros)
  expect_equal(result_dense[3, ], c(0, 0))
  
  # Check row 3 (cell_id 3, from expression row 3)
  expect_equal(result_dense[4, ], c(3.0, 0))
  
  # Check rows 4-99 are all zero
  expect_equal(sum(result_dense[5:100, ]), 0)
  
  # Check row 100 (cell_id 100, from expression row 4)
  expect_equal(result_dense[101, ], c(0, 4.0))
})

test_that("expandMatrixToCellIDs preserves sparsity", {
  # Create a sparse matrix with only a few entries
  i <- c(1, 3)
  j <- c(1, 2)
  x <- c(5.0, 6.0)
  
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(3, 2))
  
  # Cell IDs with large gaps
  all_cell_ids <- c(0, 10, 1000)
  
  result <- expandMatrixToCellIDs(mat, all_cell_ids)
  
  # Should be 1001 x 2 sparse matrix
  expect_equal(nrow(result), 1001)
  expect_equal(ncol(result), 2)
  
  # Should only have 2 non-zero entries
  expect_equal(length(result@x), 2)
  
  # Verify the values are still correct
  result_dense <- as.matrix(result)
  expect_equal(result_dense[1, 1], 5.0)  # cell_id 0, from expr row 1
  expect_equal(result_dense[1001, 2], 6.0)  # cell_id 1000, from expr row 3
})

test_that("expandMatrixToCellIDs sequential cell_ids works as identity mapping plus padding", {
  # Create a simple sparse matrix
  i <- c(1, 2)
  j <- c(1, 2)
  x <- c(1.0, 2.0)
  
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(2, 2))
  
  # Sequential cell_ids starting from 0
  all_cell_ids <- c(0, 1)
  
  result <- expandMatrixToCellIDs(mat, all_cell_ids)
  
  # Should have 2 rows (max cell_id is 1, so 0-1)
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  
  # Values should be in correct positions
  result_dense <- as.matrix(result)
  expect_equal(result_dense[1, 1], 1.0)
  expect_equal(result_dense[2, 2], 2.0)
})

