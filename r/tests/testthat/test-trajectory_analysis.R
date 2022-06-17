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


test_that("generateGraphData converts Seurat object to Monocle3 cell_data_set object", {
  data <- mock_scdata()

  cell_data <- generateGraphData(data)

  expect_s4_class(cell_data, "cell_data_set")
})

test_that("runGenerateTrajectoryGraph returns an object of class character", {
  data <- mock_scdata()

  node_umap_coords <- runGenerateTrajectoryGraph("", data)

  expect_type(node_umap_coords, "character")
})


test_that("runGenerateTrajectoryGraph json output has the expected list format when converted back to list", {
  data <- mock_scdata()

  node_umap_coords <- runGenerateTrajectoryGraph("", data)

  # convert back to list from json
  item <- RJSONIO::fromJSON(node_umap_coords)
  expect_named(item, c("nodes", "umap"))
  expect_named(item$nodes[[1]], c("x", "y", "node_id", "connected_nodes"))
  expect_named(item$umap[[1]], c("x", "y"))
  expect_type(item$nodes[[1]]$x, "double")
  expect_type(item$nodes[[1]]$y, "double")
  expect_type(item$nodes[[1]]$node_id, "character")
  expect_type(item$nodes[[1]]$connected_nodes, "character")
  expect_type(item$umap[[1]], "double")
})


