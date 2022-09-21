mock_req <- function() {
  cellSets <- list(
    children = list(
      louvain1 = list(cellIds = 0:39),
      louvain2 = list(cellIds = 40:79)
    )
  )
  req <- list(body = list(nGenes = 5, cellSets = cellSets))
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

test_that("Marker heatmap returns same sized matrices", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runMarkerHeatmap(req, data)

  expect_equal(lapply(res$rawExpression, length),
               lapply(res$truncatedExpression, length))

})


test_that("Marker Heatmap returns appropiate format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runMarkerHeatmap(req, data)

  expect_equal(names(res), c("order", "stats", "rawExpression", "truncatedExpression"))

  # number of rows in sparse matrix equals number of cells
  expect_equal(res$rawExpression$size[1], ncol(data))
  expect_equal(res$truncatedExpression$size[1], ncol(data))


  # returning only at most limit number of genes
  expect_lte(
    length(res$stats),
    req$body$nGenes * length(req$body$cellSets$children)
  )
})

test_that("Marker Heatmap nFeatures works appropiately", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 2

  res <- runMarkerHeatmap(req, data)

  # number of rows is number of cells
  expect_equal(res$rawExpression$size[1], ncol(data))

  # returning only at most limit number of genes
  expect_lte(length(res$stats), req$body$nGenes * length(req$body$cellSets$children))
})

test_that("Only one group throws error", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 5
  req$body$cellSets$children <- req$body$cellSets$children[1]

  expect_error(runMarkerHeatmap(req, data))
})
