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

  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[1]]))), max(data$cells_id) + 1)
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



test_that("Expression task continues if gene doesn't exist and returns appropriate number of results", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$genes <- append(req$body$genes, "aaa")

  res <- runExpression(req, data)

  exp <- list(
    raw = res$rawExpression,
    trunc = res$truncatedExpression,
    z = res$zScore
  )

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

  expect_equal(length(res$order), 1)
  expect_equal(length(res$stats), 1)

  expect_equal(unique(unlist(lapply(exp, function(x) x$size[[2]]))), 1)
})

test_that("If the max truncated expression value is 0, iterates correctly to find a non zero value", {
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
