runTrajectoryAnalysis <- function(req,data){
    cell_data <- SeuratWrappers::as.cell_data_set(data)
    cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
    cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

    cells_start <- req$body$cells_id
    cell_bcds <- rownames(data@meta.data)[match(cells_start,data@meta.data$cells_id)]

    cell_data <- monocle3::order_cells(cell_data,root_cells=cell_bcds)

    pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)
    pseudotime$cells_id <- data@meta.data$cells_id
    pseudotime <- pseudotime[order(pseudotime$cells_id), ]
    pseudotime <- pseudotime %>%
        tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
        select(-cells_id)
    return(unname(pseudotime))
}



