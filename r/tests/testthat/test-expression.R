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

test_that("Expression task returns appropiate number of cells and genes.", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runExpression(req, data)
  expect_equal(length(res), 2)

  expect_equal(nrow(res$rawExpression), nrow(res$truncatedExpression))
  expect_equal(ncol(res$rawExpression), ncol(res$truncatedExpression))

  expect_equal(nrow(res$rawExpression), ncol(data))
  expect_equal(ncol(res$rawExpression), 2)
})

test_that("Expression task continues if gene doesnt exist and returns appropiate number of results", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- append(req$body$genes, "aaa")

  res <- runExpression(req, data)
  expect_equal(ncol(res$rawExpression), 2)
})

test_that("Expression task works with one gene", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("MS4A1")

  res <- runExpression(req, data)
  expect_equal(ncol(res$rawExpression), 1)
  expect_equal(nrow(res$rawExpression), ncol(data))
})

test_that("runExpression throws an error if request only non existing genes", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- list("blah")

  expect_error(runExpression(req, data), "Gene\\(s\\): blah not found!")
})
