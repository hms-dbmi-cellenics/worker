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
#' @param data SeuratObject
#'
#' @return a json with nodes and umap coordinates
#' @export
runGenerateTrajectoryGraph <- function(req, data) {
  cell_data <- generateGraphData(data)
  node_coords <- t(cell_data@principal_graph_aux[["UMAP"]]$dp_mst)
  umap_coords <- as.data.frame(SingleCellExperiment::reducedDims(cell_data)[["UMAP"]])

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

  # umap data
  # fill in the NULL values for filtered cells
  umap_coords <- fillNullForFilteredCells(umap_coords, data)

  # create list
  colnames(umap_coords) <- c("x", "y")
  umap_coords_list <- lapply(asplit(umap_coords, 1), as.list)

  # node + umap data
  node_umap_coords <- list(nodes = node_coords_list, umap = umap_coords_list)

  # convert list to json
  node_umap_coords <- RJSONIO::toJSON(node_umap_coords)
  return(node_umap_coords)
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
#' earlier pseudotime growing into cells in later pseudotime, following the trajectory.
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
  root_nodes <- req$body$rootNodes

  cell_data <- generateGraphData(data)
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
generateGraphData <- function(data) {
  cell_data <- SeuratWrappers::as.cell_data_set(data)

  set.seed(ULTIMATE_SEED)

  cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
  cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

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
