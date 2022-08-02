#' Generate UMAP and node coordinates for the initial trajectory analysis plot
#'
#' Returns a json object with UMAP and node coordinates,
#' to be used by the UI for the initial plot (embedding + trajectory nodes)
#'
#' This represents the first step of the trajectory analysis.
#' It allows the creation of the initial trajectory analysis plot,
#' which is composed of an embedding and starting nodes (vertices in the lines).
#' This plot is used to select the root nodes among the starting nodes.
#' The root nodes will be then used in the following step of the trajectory
#' analysis for pseudotime calculation.
#'
#' @param req list of configuration parameters. Not used here.
#' @param data SeuratObject
#'
#' @return a json with nodes and umap coordinates
#' @export
runGenerateTrajectoryGraph <- function(req, data) {
  cell_data <- generateGraphData(
    req$body$embedding,
    req$body$embedding_settings,
    req$body$clustering_settings,
    data
  )

  node_coords <- t(cell_data@principal_graph_aux[["UMAP"]]$dp_mst)

  # node coordinates
  # get connected nodes
  connected_nodes <- list()
  for (node in rownames(node_coords)) {
    node_id <- which(rownames(node_coords) == node)
    connected_nodes_obj <- cell_data@principal_graph[["UMAP"]][[node_id]][[1]]
    connected_nodes[[node]] <- as.list(names(connected_nodes_obj))
  }

  # create list
  colnames(node_coords) <- c("x", "y")
  node_coords_list <- lapply(asplit(node_coords, 1), as.list)
  for (i in 1:length(node_coords_list)) {
    node_coords_list[[i]]["node_id"] <- names(node_coords_list)[i]
    node_coords_list[[i]]["connected_nodes"] <- list(connected_nodes[[i]])
  }

  # node + umap data
  root_nodes <- list(nodes = node_coords_list)
  return(root_nodes)
}


#' Calculate pseudotime
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
#' @param req {body: {
#'               rootNodes: root nodes ids. Determines the root nodes of the trajectory
#'              }
#'            }
#' @param data SeuratObject
#'
#' @return a tibble with pseudotime values
#' @export
runTrajectoryAnalysis <- function(req, data) {
  cell_data <- generateGraphData(
    req$body$embedding,
    req$body$embedding_settings,
    req$body$clustering_settings,
    data
  )

  root_nodes <- req$body$rootNodes

  cell_data <- monocle3::order_cells(cell_data, reduction_method = "UMAP", root_pr_nodes = root_nodes)

  pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)

  # fill in the NULL values for filtered cells
  pseudotime <- fillNullForFilteredCells(pseudotime, data)

  return(unname(pseudotime))
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
generateGraphData <- function(embedding_data, embedding_settings, clustering_settings, data) {
  set.seed(ULTIMATE_SEED)

  Seurat::DefaultAssay(data) <- "RNA"

  clustering_method <- clustering_settings$method
  embedding_method <- embedding_settings$method

  # Clustering resolution can only be used by monocle if the clustering method is leiden
  clustering_resolution <- NULL
  if(clustering_method == "leiden") {
    clustering_resolution <- clustering_settings$methodSettings[[clustering_method]]$resolution
  }

  data <- assignEmbedding(embedding_data, data)

  surat_to_monocle_method_map <- list(
    umap = "UMAP",
    tsne = "tSNE"
  )

  cell_data <- SeuratWrappers::as.cell_data_set(data)

  message("Calculating trajectory graph...")
  cell_data <- monocle3::cluster_cells(
    cds = cell_data,
    cluster_method = clustering_method,
    reduction_method = surat_to_monocle_method_map[[embedding_method]],
    resolution = clustering_resolution
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
