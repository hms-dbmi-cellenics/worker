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


test_that("Marker heatmap returns same sized matrices", {
  data <- mock_scdata()
  req <- mock_req()

  # dims of returned expression data includes all the downsampled
  expected_sizes <-
    c(
      max(req$body$cellIds) + 1,
      req$body$nGenes * length(req$body$cellSets$children)
    )

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  sizes <- list(
    res$rawExpression$size,
    res$truncatedExpression$size,
    res$zScore$size
  )

  sizes <- lapply(sizes, unlist)

  expect_true(all(unlist(
    lapply(sizes, all.equal, expected_sizes)
  )))
})


test_that("Marker Heatmap returns appropiate format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  expect_equal(
    names(res),
    c("orderedGeneNames", "stats", "rawExpression", "truncatedExpression", "zScore")
  )

  # number of rows in sparse matrix equals number of cells
  expect_equal(
    unlist(res$rawExpression$size)[1],
    max(req$body$cellIds) + 1
  )
  expect_equal(
    unlist(res$truncatedExpression$size)[1],
    max(req$body$cellIds) + 1
  )


  # returning only at most limit number of genes
  expect_lte(
    length(res$stats),
    req$body$nGenes * length(req$body$cellIds)
  )
})


test_that("Marker Heatmap nFeatures works appropiately", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 2

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  # number of rows is number of cells
  expect_equal(
    unlist(res$rawExpression$size[1]),
    max(req$body$cellIds) + 1
  )

  # returning only at most limit number of genes per cellset
  expect_true(all(
    unlist(lapply(res$stats, length)) <= req$body$nGenes * length(req$body$cellSets$children)
  ))
})


test_that("all stats contain the same number of elements", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runMarkerHeatmap(req, data)
  withr::defer(cleanupMarkersCache())

  expect_true(all(lapply(res$stats, length) == length(res$stats[[1]])))
})


test_that("Only one group throws error", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$nGenes <- 5
  req$body$cellSets$children <- req$body$cellSets$children[1]

  expect_error(runMarkerHeatmap(req, data))
  withr::defer(cleanupMarkersCache())
})
