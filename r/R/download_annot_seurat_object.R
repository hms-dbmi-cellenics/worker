#' Annotate Seurat object
#'
#' This function adds information about embeddings and clustering to the processed Seurat object.
#' Embedding are obtained and added to the Seurat object using the assignEmbedding function.
#' Cellsets object is received from the Python worker and parsed to extract
#' relevant information about louvain/leiden clusters, custom cellsets,
#' and cluster annotations. This information is added to the Seurat object
#' to provide a more comprehensive representation of the data.
#'
#'
#' @param req Request object that contains the cellsets object and embedding data.
#' @param data Seurat object
#'
#' @return Path to the saved Seurat object.
#' @export
#'
DownloadAnnotSeuratObject <- function(req, data) {
  cell_sets <- req$body$cellSets

  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  if (!req$body$isObj2s) {
    embedding_data <- req$body$embedding
    embedding_method <- req$body$embeddingMethod
    data <- assignEmbedding(embedding_data, data, embedding_method)
  }

  saveRDS(data, RDS_PATH)
  return(RDS_PATH)
}
