mock_req <- function() {
  cellSets <- list(
    children = list(
      louvain1 = list(cellIds = c(0:30, 40:65)),
      louvain2 = list(cellIds = c(66:70, 100:117))
    )
  )
  cellIds <- c(0:20, 40:55)

  req <- list(body = list(nGenes = 5, cellSets = cellSets, cellIds = cellIds))
  return(req)
}


mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- c(0:30, 40:70, 100:117)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )
  return(pbmc_small)
}


test_that("Marker heatmap returns orderedGeneNames", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  expect_equal(
    names(res),
    c("orderedGeneNames")
  )

  # should have genes for each cellset
  expect_equal(
    length(res$orderedGeneNames),
    req$body$nGenes * length(req$body$cellSets$children)
  )
})


test_that("Marker Heatmap nFeatures works appropriately", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 2

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  # returning only nGenes genes per cellset
  expect_equal(
    length(res$orderedGeneNames),
    req$body$nGenes * length(req$body$cellSets$children)
  )
})


test_that("Only one group throws error", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 5
  req$body$cellSets$children <- req$body$cellSets$children[1]

  expect_error(runMarkerHeatmap(req, data))
  withr::defer(cleanupMarkersCache())
})
