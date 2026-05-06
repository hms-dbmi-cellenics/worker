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


mock_scdata_marker_heatmap <- function(use_bpcells = FALSE) {
  scdata <- mock_scdata(use_bpcells = use_bpcells)
  scdata$cells_id <- c(0:30, 40:70, 100:117)
  return(scdata)
}

test_that("Marker heatmap works with bpcells", {
  data <- mock_scdata_marker_heatmap(use_bpcells = TRUE)
  req <- mock_req()

  expect_no_error(runMarkerHeatmap(req, data))
})

test_that("Marker heatmap returns orderedGeneNames", {
  data <- mock_scdata_marker_heatmap()
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
  data <- mock_scdata_marker_heatmap()
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
  data <- mock_scdata_marker_heatmap()
  req <- mock_req()
  req$body$nGenes <- 5
  req$body$cellSets$children <- req$body$cellSets$children[1]

  expect_error(runMarkerHeatmap(req, data))
  withr::defer(cleanupMarkersCache())
})
