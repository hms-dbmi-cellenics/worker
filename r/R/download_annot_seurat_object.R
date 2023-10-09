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
  embedding_data <- req$body$embedding

  data <- add_cellsets(data, cell_sets)
  data <- assignEmbedding(embedding_data, data)

  saveRDS(data, RDS_PATH)

  return(RDS_PATH)
}



add_cellsets <- function(scdata, cellsets) {

  for (i in seq_along(cellsets)) {

    cellset <- cellsets[[i]]

    # get children
    children <- cellset$children
    if (!length(children)) next()

  if (uuid::UUIDvalidate(cellset$key)) {
        key <- cellset$name
    } else {
        key <- cellset$key
    }

    scdata[[key]] <- NA_character_



    for (j in seq_along(children)) {
      child <- children[[j]]

      # use name (aka user supplied label) as value
      value <- child$name

      # get associated 'cells_id'
      cells_id <- unlist(child$cellIds)

      # add value for relevant 'cells_id'
      is.cells_id <- scdata$cells_id %in% cells_id
      scdata@meta.data[is.cells_id, key] <- value
    }
  }
  return(scdata)
}
