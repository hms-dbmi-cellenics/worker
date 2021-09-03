runTrajectoryAnalysis <- function(req, data){

    # Available reduction method: "UMAP", "tSNE", "PCA", "LSI", "Aligned"
    # Because we use a different case, we have to translate them
    # method <- switch(
    #     req$body$method,
    #     "umap" = "UMAP",
    #     "tsne" = "tSNE"
    # )

    cell_data <- SeuratWrappers::as.cell_data_set(data)
    cell_data <- monocle3::cluster_cells(cds = cell_data, reduction_method = "UMAP")
    cell_data <- monocle3::learn_graph(cell_data, use_partition = TRUE)

    cells_start <- req$body$cell_ids
    cell_bcds <- rownames(data@meta.data)[match(cells_start, data@meta.data$cells_id)]

    cell_data <- monocle3::order_cells(cell_data,root_cells=cell_bcds)

    pseudotime <- as.data.frame(cell_data@principal_graph_aux@listData$UMAP$pseudotime)
    pseudotime$cells_id <- data@meta.data$cells_id
    pseudotime <- pseudotime[order(pseudotime$cells_id), ]
    pseudotime <- pseudotime %>%
        tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
        select(-cells_id)

    # Disabled until path algorithm can be confirmed
    # ica_space_df <- t(cell_data@principal_graph_aux$UMAP$dp_mst) %>%
    #   as.data.frame() %>%
    #   dplyr::select(x = 1, y = 2) %>%
    #   dplyr::mutate(cluster_name = rownames(.))

    # dp_mst <- cell_data@principal_graph$UMAP

    # edge_df <- dp_mst %>%
    #   igraph::as_data_frame() %>%
    #   dplyr::select(source = "from", target = "to") %>%
    #   dplyr::left_join(ica_space_df %>%
    #                      dplyr::select(source="cluster_name",
    #                                     source_x="x",
    #                                     source_y="y"),
    #                    by = "source") %>%
    #   dplyr::left_join(ica_space_df %>%
    #                      dplyr::select(target="cluster_name",
    #                                     target_x="x",
    #                                     target_y="y"),
    #                    by = "target")

    # edges <- list()
    # for (i in 1:nrow(edge_df)) {
    #       edges <- append(edges,
    #         list(
    #             list(
    #                 x1 = edge_df[i, "source_x"],
    #                 y1 = edge_df[i, "source_y"],
    #                 x2 = edge_df[i, "target_x"],
    #                 y2 = edge_df[i, "target_y"]
    #             )
    #         )
    #     )
    # }

    return(
        list(
            pseudotime = unname(pseudotime),
            # graph = edges -- disabled until path algorithm can be confirmed 
            graph = list()
        )
    )
}