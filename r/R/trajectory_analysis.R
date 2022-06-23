runTrajectoryAnalysis <- function(req, data) {
  root_nodes <- req$body$rootNodes

  cell_data <- generateGraphData(data)
  cell_data <- monocle3::order_cells(cell_data, reduction_method = "UMAP", root_pr_nodes = root_nodes)

  pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)
  pseudotime$cells_id <- data@meta.data$cells_id
  pseudotime <- pseudotime[order(pseudotime$cells_id), ]
  pseudotime <- pseudotime %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    dplyr::select(-cells_id)
  return(unname(pseudotime))
}

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
  umap_coords$cells_id <- data@meta.data$cells_id
  umap_coords <- umap_coords[order(umap_coords$cells_id), ]
  umap_coords <- umap_coords %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    dplyr::select(-cells_id)
  # create list
  colnames(umap_coords) <- c("x", "y")
  umap_coords_list <- lapply(asplit(umap_coords, 1), as.list)

  # node + umap data
  node_umap_coords <- list(nodes = node_coords_list, umap = umap_coords_list)

  # convert list to json
  node_umap_coords <- RJSONIO::toJSON(node_umap_coords)
  return(node_umap_coords)
}


generateGraphData <- function(data) {
  cell_data <- SeuratWrappers::as.cell_data_set(data)

  set.seed(42)

  cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
  cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

  return(cell_data)
}
