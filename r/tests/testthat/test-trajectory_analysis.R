mock_scdata <- function(filt_cell_id = "") {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  if (all(filt_cell_id != "")) {
    keep_cells_id <- which(!pbmc_small$cells_id %in% filt_cell_id)
    keep_cells <- names(pbmc_small$cells_id[keep_cells_id])
    pbmc_small <- subset(pbmc_small, cells = keep_cells)
  }
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )

  # scale and PCA
  pbmc_small <- Seurat::NormalizeData(pbmc_small, normalization.method = "LogNormalize", verbose = FALSE)
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::RunPCA(pbmc_small, verbose = FALSE, npcs = 10)
  pbmc_small@misc[["active.reduction"]] <- "pca"

  # run UMAP
  npcs <- get_npcs(pbmc_small)
  pbmc_small <- suppressWarnings(Seurat::RunUMAP(pbmc_small, dims = 1:npcs, verbose = FALSE))

  return(pbmc_small)
}

get_explained_variance <- function(scdata) {
  # Compute explained variance for plotting and numPCs estimation.
  # It can be computed from pca or other reductions such as mnn
  if (scdata@misc[["active.reduction"]] == "mnn") {
    var_explained <- scdata@tools$`SeuratWrappers::RunFastMNN`$pca.info$var.explained
  } else {
    eig_values <- (scdata@reductions$pca@stdev)^2
    var_explained <- eig_values / sum(eig_values)
  }
  return(var_explained)
}

get_npcs <- function(scdata, var_threshold = 0.85, max_npcs = 30) {
  # estimates the number of PCs to use in data integration and embeddings,
  # using accumulated explained variance
  var_explained <- get_explained_variance(scdata)
  npcs <- min(which(cumsum(var_explained) >= var_threshold))
  return(min(npcs, max_npcs, na.rm = TRUE))
}


mock_starting_nodes_req <- function(data) {
  mock_embedding_data <- get_mock_embedding_data(data)
  mock_embedding_settings <- get_mock_embedding_settings()
  mock_clustering_settings <- get_mock_clustering_settings()

  result <- list(
    body = list(
      embedding = mock_embedding_data,
      embedding_settings = mock_embedding_settings,
      clustering_settings = mock_clustering_settings
    )
  )
}

mock_pseudotime_req <- function(data) {
  mock_embedding_data <- get_mock_embedding_data(data)
  mock_embedding_settings <- get_mock_embedding_settings()
  mock_clustering_settings <- get_mock_clustering_settings()

  result <- list(
    body = list(
      embedding = mock_embedding_data,
      embedding_settings = mock_embedding_settings,
      clustering_settings = mock_clustering_settings,
      root_nodes = c(0, 1, 2)
    )
  )

}

get_mock_embedding_data <- function(mock_data) {
  embedding_data <- mock_data@reductions$umap@cell.embeddings
  embedding_data <- unname(embedding_data)
  embedding_data <- lapply(
    seq_len(nrow(embedding_data)),
    function(i) embedding_data[i,]
  )
}

get_mock_embedding_settings <- function() {
  result <- list (
    Etag = "mockEmbeddingETag",
    method = "umap",
    methodSettings = list (
      distanceMetric = "cosine",
      minimumDistance = 0.3
    )
  )
}

get_mock_clustering_settings <- function() {
  result <- list (
    method = "louvain",
    resolution = 0.8
  )
}

test_that("generateTrajectoryGraph converts Seurat object to Monocle3 cell_data_set object", {
  data <- mock_scdata()

  mock_embedding_data <- get_mock_embedding_data(data)
  mock_embedding_settings <- get_mock_embedding_settings()
  mock_clustering_settings <- get_mock_clustering_settings()

  cell_data <- suppressWarnings(
    generateTrajectoryGraph(
      mock_embedding_data,
      mock_embedding_settings,
      mock_clustering_settings,
      data
    )
  )

  expect_s4_class(cell_data, "cell_data_set")
})


test_that("runTrajectoryAnalysisStartingNodesTask output has the expected format", {
  data <- mock_scdata()
  req <- mock_starting_nodes_req(data)

  root_nodes <- suppressWarnings(runTrajectoryAnalysisStartingNodesTask(req, data))

  expect_named(root_nodes, c("connectedNodes", "x", "y"))
  expect_type(root_nodes$x[[1]], "double")
  expect_type(root_nodes$y[[1]], "double")
  expect_type(root_nodes$connectedNodes[[1]], "double")

  # Second element is a list because it has length 1, so this way we ensure that
  # it is encoded in json as a list too
  expect_type(root_nodes$connectedNodes[[2]], "list")
})


test_that("runTrajectoryAnalysisStartingNodesTask outputs the correct number of nodes", {
  data <- mock_scdata()
  mock_embedding_data <- get_mock_embedding_data(data)
  mock_embedding_settings <- get_mock_embedding_settings()
  mock_clustering_settings <- get_mock_clustering_settings()

  req <- mock_starting_nodes_req(data)

  cell_data <- suppressWarnings(
    generateTrajectoryGraph(
      mock_embedding_data,
      mock_embedding_settings,
      mock_clustering_settings,
      data
    )
  )

  root_nodes <- suppressWarnings(runTrajectoryAnalysisStartingNodesTask(req, data))
  expected_length <- nrow(t(cell_data@principal_graph_aux[["UMAP"]]$dp_mst))

  expect_equal(expected_length, length(root_nodes$connectedNodes))
  expect_equal(expected_length, length(root_nodes$x))
  expect_equal(expected_length, length(root_nodes$y))
})


test_that("runTrajectoryAnalysisPseudoTimeTask works", {
  data <- mock_scdata()
  req <- mock_pseudotime_req(data)

  result <- suppressWarnings(runTrajectoryAnalysisPseudoTimeTask(req, data))

  expect_equal(length(data$cells_id), length(result$pseudotime))
  expect_true(is.numeric(result$pseudotime))
})

test_that('runTrajectoryAnalysisPseudoTimeTask fails if root_node is empty', {
  data <- mock_scdata()
  req <- mock_pseudotime_req(data)

  # Set root nodes to empty vector
  req$body$root_nodes <- c()

  expect_error(runTrajectoryAnalysisPseudoTimeTask(req, data), "No root nodes were selected for the analysis.")
})
