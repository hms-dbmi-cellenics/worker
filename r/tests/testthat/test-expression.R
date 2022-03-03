mock_req <- function() {
  req <- list(body = list(genes = list("MS4A1", "CD79B")))
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

test_that("Expression task returns appropiate number and names of genes.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)

  expect_equal(length(res), 2)
  expect_equal(as.list(names(res)), req$body$genes)
})

test_that("Expression task returns appropiate number of cells.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)
  raw_expression <- res[[1]]$rawExpression$expression
  truncated_expression <- res[[1]]$truncatedExpression$expression


  expect_equal(length(raw_expression), length(truncated_expression))

  expect_equal(length(raw_expression), ncol(data))
})

test_that("Expression result is properly formated.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)
  raw_expression_result <- res[[1]]$rawExpression
  truncated_expression_result <- res[[1]]$truncatedExpression
  expect_true(!all(is.na(raw_expression_result$expression)))

  expect_equal(raw_expression_result$mean, mean(raw_expression_result$expression, na.rm = TRUE))
  expect_equal(raw_expression_result$stdev, sd(raw_expression_result$expression, na.rm = TRUE))
  expect_equal(truncated_expression_result$min, min(truncated_expression_result$expression, na.rm = TRUE))
  expect_equal(truncated_expression_result$max, max(truncated_expression_result$expression, na.rm = TRUE))
})



test_that("Expression task continues if gene doesnt exist and returns appropiate number of results", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- append(req$body$genes, "aaa")

  res <- runExpression(req, data)
  expect_equal(length(res), 2)
})

test_that("Expression task works with one gene", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("MS4A1")

  res <- runExpression(req, data)
  expect_equal(length(res), 1)
  expect_equal(length(res[[1]]$rawExpression$expression), ncol(data))
})

test_that("If the max truncated expression value is 0, iterates correctly to find a non zero value", {
  data <- mock_scdata()
  req <- mock_req()

  data@assays$RNA@data["MS4A1",] <- 0
  data@assays$RNA@data["MS4A1",1] <- 5

  res <- runExpression(req, data)
  expect_false(all(res$MS4A1$truncatedExpression$expression==0) )
})

test_that("runExpression throws an error if request only non existing genes", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("blah")

  expect_error(runExpression(req, data), "Gene\\(s\\): blah not found!")
})
