#' Returns expression values for selected genes
#'
#' @param genes - Must have names and input(ensmbl ids)
#' @param data
#'
#' @return
#' @export
#'
#' @examples
getExpressionValues <- function(genes,data){
    quantile_threshold <- 0.95
    #
    #Get the expression values for those genes in the corresponding matrix.
    geneExpression <-data@assays$RNA@data[unique(genes$input),,drop=FALSE]
    geneExpression <- as.data.frame(t(as.matrix(geneExpression)))
    geneExpression$cells_id <- data@meta.data$cells_id
    geneExpression <- geneExpression[ order(geneExpression$cells_id), ]
    geneExpression <- geneExpression %>%
        tidyr::complete(cells_id = seq(0,max(data@meta.data$cells_id))) %>%
        select(-cells_id)
    # worried about duplicate gene row.names in @data
    symbol_idx <- match(colnames(geneExpression), genes$input)
    colnames(geneExpression) <- genes$name[symbol_idx]
    adjGeneExpression <- as.data.frame(apply(geneExpression,2,FUN=function(x){
        lim <- as.numeric(quantile(x,quantile_threshold,na.rm=TRUE))
        i <- 0.01
        while (lim==0 & i+quantile_threshold<=1) {
            lim <- as.numeric(quantile(x,quantile_threshold + i,na.rm=TRUE))
            i<-i+0.01
        }
        print(i)
        return(pmin(x,lim))
    }))
    return(list(rawExpression = geneExpression,truncatedExpression = adjGeneExpression))
}

#' Compute clusters and return object with clusters
#'
#' @param algorithm
#' @param resolution
#' @param data
#'
#' @return
#'
#' @examples
getClusters <- function(type, resolution, data) {
    algorithm <- list("louvain" = 1, "leiden" = 4)[[type]]
    # To run clustering, we need to identify the active.reduction that is used in FindNeighbors.
    if ("active.reduction" %in% names(data@misc)) {
        active.reduction <- data@misc[["active.reduction"]]
    } else {
        active.reduction <- "pca"
    }

    if (type == "leiden") {
        # emulate FindClusters, which overwrites seurat_clusters slot and meta.data column
        g <- getSNNiGraph(data)
        clus_res <- igraph::cluster_leiden(g, "modularity", resolution_parameter = resolution)
        clusters <- clus_res$membership
        names(clusters) <- clus_res$names
        clusters <- clusters[colnames(data)]
        data$seurat_clusters <- data@meta.data[, res_col] <- factor(clusters - 1)
    } else {
        graph.name <- paste0(DefaultAssay(data), "_snn")
        if (!graph.name %in% names(data)) {
            data <- Seurat::FindNeighbors(data, k.param = 20, annoy.metric = "cosine", verbose = FALSE, reduction = active.reduction)
        }
        data <- FindClusters(data, resolution = resolution, verbose = FALSE, algorithm = algorithm)
    }
    print(resolution)
    return(data)
}
