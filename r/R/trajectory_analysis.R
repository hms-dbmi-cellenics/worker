runTrajectoryAnalysis <- function(req,data){
    root_nodes <- req$body$rootNodes

    cell_data <- generateGraphData(data)
    cell_data <- monocle3::order_cells(cell_data,reduction_method="UMAP", root_pr_nodes = root_nodes)

    #monocle3::plot_cells(cell_data,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)

    pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)
    pseudotime$cells_id <- data@meta.data$cells_id
    pseudotime <- pseudotime[order(pseudotime$cells_id), ]
    pseudotime <- pseudotime %>%
        tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
        select(-cells_id)
    return(unname(pseudotime))
}

runGenerateTrajectoryGraph <- function(req,data){
    cell_data <- generateGraphData(data)
    node_coords <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst)
    umap_coords <- as.data.frame(SingleCellExperiment::reducedDims(cds)[["UMAP"]])

    #TO DO
    #return the node and umap coords according to the design doc
    #for the umap coords, fill in the NULL values for filtered cells
    #Response format:
    #   {
    #    nodes:  {
    #        node_id: {
    #            x,
    #            y,
    #            node_id,
    #            connected_nodes = [node_id, node_id]
    #        }
    #    },
    #    Umap: [{x, y}, {x, y, {x, y}],
    #    }
    #
}

generateGraphData <- function(data){
    cell_data <- SeuratWrappers::as.cell_data_set(data)

    set.seed(42)

    cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
    cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

    return(cell_data)
}
