  # Currently Monocle3 only supports trajectory analysis for UMAP embeddings,
  # that is why there is only 1 method in the mapping below. The embedding method
  # that is used for embedding and trajectory generation is hard coded in the UI.
  SEURAT_TO_MONOCLE_METHOD_MAP <- list(
    umap = "UMAP"
  )

#' Generate node coordinates for the initial trajectory analysis plot
#'
#' Returns a list containing node coordinates to the python worker
#'
#' This represents the first step of the trajectory analysis.
#' It allows the creation of the initial trajectory analysis plot,
#' which is composed of an embedding and starting nodes (vertices in the lines).
#' This plot is used to select the root nodes among the starting nodes.
#' The root nodes will be then used in the following step of the trajectory
#' analysis for pseudotime calculation.
#'
#' @param req {
#'            body: {
#'               embedding: numeric vector containing a pair of embedding coordinates (e.g. num [1:2] -1, -2)
#'               embedding_settings: {
#'                  method: Embedding method (.e.g umap)
#'                  methodSettings: An object containing settings specific for each embedding type
#'               },
#'               clustering_settings: {
#'                  method: Clustering method (e.g. louvain),
#'                  resolution: Clustering resolution
#'               },
#'              }
#'            }
#' @param data SeuratObject
#'
#' @return a list containing nodes coordinates, connected nodes and the node_id
#' @export
runTrajectoryAnalysisStartingNodesTask <- function(req, data) {
  cell_data <- generateTrajectoryGraph(
    req$body$embedding,
    req$body$embedding_settings,
    req$body$clustering_settings,
    req$body$cell_ids,
    data
  )

  seurat_embedding_method <- req$body$embedding_settings$method
  monocle_embedding_method <- SEURAT_TO_MONOCLE_METHOD_MAP[[seurat_embedding_method]]

  node_coords <- t(cell_data@principal_graph_aux[[monocle_embedding_method]]$dp_mst)

  # node coordinates
  # get connected nodes
  connected_nodes <- list()

  for (node_id in seq_along(rownames(node_coords))) {

    # Get connected node index vector
    current_connected_nodes <- cell_data@principal_graph[[monocle_embedding_method]][[node_id]][[1]]

    # Keep only those that are higher than current node
    current_connected_nodes <- as.integer(current_connected_nodes[current_connected_nodes > node_id])

    # Shift by 1 to use 0-based indexes
    current_connected_nodes <- current_connected_nodes - 1

    connected_nodes[[node_id]] <- ensure_is_list_in_json(current_connected_nodes)
  }

  return(
    list(
      connectedNodes = connected_nodes,
      x = unname(node_coords[, 1]),
      y = unname(node_coords[, 2])
    )
  )
}


#' Calculate trajectory analysis pseudotime
#'
#' Order the cells and generate an array of pseudotime values, based on
#' the node ids of the root nodes.
#'
#' Pseudotime is a value that shows the relative progression
#' of a cell through a path in gene expression space.
#' These values will be used to plot the cells along a continuous path
#' that represents the evolution of the process.
#' This plot represents the final plot of the trajectory analysis, which shows
#' a map of how cells expression changes, starting from cells with smaller
#' pseudotime values to cells with larger pseudotimes, along a trajectory.
#'
#' @param req {
#'            body: {
#'               embedding: numeric vector containing a pair of embedding coordinates (e.g. num [1:2] -1, -2)
#'               embedding_settings: {
#'                  method: Embedding method (.e.g umap)
#'                  methodSettings: An object containing settings specific for each embedding type
#'               },
#'               clustering_settings: {
#'                  method: Clustering method (e.g. louvain),
#'                  resolution: Clustering resolution
#'               },
#'               root_nodes: root nodes ids. Determines the root nodes of the trajectory
#'              }
#'            }
#' @param data SeuratObject
#'
#' @return a tibble with pseudotime values
#' @export
runTrajectoryAnalysisPseudoTimeTask <- function(req, data) {
  if(length(req$body$root_nodes) == 0) {
    stop(
      generateErrorMessage(
        error_codes$EMPTY_ROOT_NODES,
        "No root nodes were selected for the analysis."
      )
    )
  }

  cell_data <- generateTrajectoryGraph(
    req$body$embedding,
    req$body$embedding_settings,
    req$body$clustering_settings,
    req$body$cell_ids,
    data
  )

  seurat_embedding_method <- req$body$embedding_settings$method
  monocle_embedding_method <- SEURAT_TO_MONOCLE_METHOD_MAP[[seurat_embedding_method]]

  node_ids <- colnames(cell_data@principal_graph_aux[[monocle_embedding_method]]$dp_mst)

  # Add 1 to indexes so they are 1-based
  root_indexes <- req$body$root_nodes + 1

  # Translate the indexes to their ids
  root_ids <- node_ids[root_indexes]

  cell_data <- monocle3::order_cells(cell_data, reduction_method = monocle_embedding_method, root_pr_nodes = root_ids)

  pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)

  # fill in the NULL values for filtered cells
  subset_data <- subsetIds(data, req$body$cell_ids)
  pseudotime <- fillNullForFilteredCells(pseudotime, subset_data)
  result <- list(pseudotime = pseudotime[[1]])
  return(result)
}


#' Convert Seurat object to Monocle3 cell_data_set object, cluster cells and learn graph
#'
#' In order to perform the trajectory analysis with Monocle3,
#' the Seurat object needs to be converted to a Monocle3 cell_data_set object.
#' After conversion, this function also learns the trajectory graph.
#'
#' @param data Seurat object
#'
#' @return a cell_data_set object with cluster and graph information stored internally
#' @export
generateTrajectoryGraph <- function(
  embedding_data,
  embedding_settings,
  clustering_settings,
  cell_ids,
  data
) {

  set.seed(ULTIMATE_SEED)

  Seurat::DefaultAssay(data) <- "RNA"

  clustering_method <- clustering_settings$method
  embedding_method <- embedding_settings$method

  # Clustering resolution can only be used by monocle if the clustering method is leiden
  clustering_resolution <- NULL
  if (clustering_method == "leiden") {
    clustering_resolution <- clustering_settings$resolution
  }

  clustering_controls <- list()
  if(embedding_method == "umap") {
    clustering_controls <- list(
      metric = embedding_settings$method_settings$distanceMetric
    )
  }

  data <- subsetIds(data, cell_ids)
  data <- assignEmbedding(embedding_data, data)

  cell_data <- SeuratWrappers::as.cell_data_set(data)

  message("Calculating trajectory graph...")
  cell_data <- monocle3::cluster_cells(
    cds = cell_data,
    cluster_method = clustering_method,
    reduction_method = SEURAT_TO_MONOCLE_METHOD_MAP[[embedding_method]],
    resolution = clustering_resolution,
    nn_control = clustering_controls,
    k = 25
  )

  cell_data <- monocle3::learn_graph(cell_data)

  return(cell_data)
}


#' Fill in the NULL values for filtered cells
#'
#' When a dataframe with cell barcodes and some associated information is generated
#' from a filtered Seurat object, it will not contain the cell barcodes that were filtered out.
#' However, the UI needs the whole array of ordered and unfiltered cell ids.
#' For this reason, there is a need to add to the dataframe the corresponding cell ids
#' for cell barcodes that were filtered out and fill the values with NULL.
#'
#' @param df data.frame with barcodes as rownames
#' @param data Seurat object with cells_id
#'
#' @return a tibble filled with NULL values for missing cell ids
#' @export
fillNullForFilteredCells <- function(df, data) {
  df$cells_id <- data@meta.data$cells_id
  max_value <- max(data@meta.data$cells_id)
  df <- df[order(df$cells_id), ]

  df <- df %>%
    tidyr::complete(cells_id = seq(0, max_value)) %>%
    dplyr::select(-cells_id)

  return(df)
}
