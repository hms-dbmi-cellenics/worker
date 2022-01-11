#
# getClusters
# returns the clusters in the shape of a dataframe with a clusters column,
# cell ids column and cell barcode as rownames.
#
# req$body has:
# type: can be "louvain"/"leiden"
# config:{
#          resolution: integer, range: 0 - 2
#         }
#
#
# We currently CANT support leiden, we need to discuss this in bioinformatics,
#  the algorithm is not working.
#
#' @export
#'
runClusters <- function(req, data) {
  resol <- req$body$config$resolution
  type <- req$body$type

  data <- getClusters(type, resol, data)
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resol))
  # In the meta data slot the clustering is stored with the resolution
  # used to calculate it
  # RNA_snn_res.#resolution
  df <- data.frame(
    cluster = data@meta.data[, res_col],
    cell_ids = data@meta.data$cells_id
  )
  # get the cell barcodes as rownames
  rownames(df) <- rownames(data@meta.data)
  return(df)
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

  # need the reduction that is used in FindNeighbors.
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

    graph.name <- paste0(Seurat::DefaultAssay(data), "_snn")
    if (!graph.name %in% names(data)) {
      data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = active.reduction)
    }
    data <- Seurat::FindClusters(data, resolution = resolution, verbose = FALSE, algorithm = algorithm)
  }

  return(data)
}

