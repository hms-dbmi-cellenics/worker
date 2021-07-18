#' Returns expression values for selected genes
#'
#' @param genes - Must have names and input(ensmbl ids)
#' @param data
#'
#' @return
#' @export
#'
#' @examples
getExpressionValues <- function(genes, data) {
  quantile_threshold <- 0.95
  #
  # Get the expression values for those genes in the corresponding matrix.
  geneExpression <- data@assays$RNA@data[unique(genes$input), , drop = FALSE]
  geneExpression <- as.data.frame(t(as.matrix(geneExpression)))
  geneExpression$cells_id <- data@meta.data$cells_id
  geneExpression <- geneExpression[order(geneExpression$cells_id), ]
  geneExpression <- geneExpression %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    select(-cells_id)
  # worried about duplicate gene row.names in @data
  symbol_idx <- match(colnames(geneExpression), genes$input)
  colnames(geneExpression) <- genes$name[symbol_idx]
  adjGeneExpression <- as.data.frame(apply(geneExpression, 2, FUN = function(x) {
    lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
    i <- 0.01
    while (lim == 0 & i + quantile_threshold <= 1) {
      lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
      i <- i + 0.01
    }
    return(pmin(x, lim))
  }))
  return(list(rawExpression = geneExpression, truncatedExpression = adjGeneExpression))
}

#' Get and Convert SNN Graph object into igraph object
#'
#' This is used to facilitate leiden clustering.
#'
#' @param data \code{Seurat} object
#'
#' @return boolean indicating if SNN Graph object exists
#'
getSNNiGraph <- function(data) {

  # check to see if we already have Seurat SNN Graph object
  snn_name <- paste0(data@active.assay, "_snn")

  # if doesn't exist, run SNN
  if (!snn_name %in% names(data)) data <- Seurat::FindNeighbors(data)

  # convert Seurat Graph object to igraph
  # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
  adj_matrix <- Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE
  )


  return(g)
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
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
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
  return(data)
}
