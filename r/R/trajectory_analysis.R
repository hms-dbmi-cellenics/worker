runGenerateTrajectoryGraph <- function(data) {
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
  umap_coords <- fill_null_for_filtered_cells(umap_coords, data)

  # create list
  colnames(umap_coords) <- c("x", "y")
  umap_coords_list <- lapply(asplit(umap_coords, 1), as.list)

  # node + umap data
  node_umap_coords <- list(nodes = node_coords_list, umap = umap_coords_list)

  # convert list to json
  node_umap_coords <- RJSONIO::toJSON(node_umap_coords)
  return(node_umap_coords)
}


runTrajectoryAnalysis <- function(req, data) {
  root_nodes <- req$body$rootNodes

  cell_data <- generateGraphData(data)
  cell_data <- monocle3::order_cells(cell_data, reduction_method = "UMAP", root_pr_nodes = root_nodes)

  pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)

  # fill in the NULL values for filtered cells
  pseudotime <- fill_null_for_filtered_cells(pseudotime, data)

  return(unname(pseudotime))
}



generateGraphData <- function(data) {
  cell_data <- SeuratWrappers::as.cell_data_set(data)

  set.seed(ULTIMATE_SEED)

  cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
  cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

  return(cell_data)
}



fill_null_for_filtered_cells <- function(df, data) {
  # add cells_id column
  df$cells_id <- data@meta.data$cells_id
  # get max value for cells_id
  max_value <- max(data@meta.data$cells_id)
  # order by cells_id
  df <- df[order(df$cells_id), ]
  # add NULL values for filtered cells and remove cells_id column
  df <- df %>%
    tidyr::complete(cells_id = seq(0, max_value)) %>%
    dplyr::select(-cells_id)

  return(df)
}
