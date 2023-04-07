DownloadAnnotSeuratObject <- function(req, data) {

  cell_sets <- req$body$cellSets
  embedding_data <- req$body$embedding

  parsed_cellsets <- parse_cellsets(cell_sets)
  data <- add_clusters(data, parsed_cellsets)

  data <- assignEmbedding(embedding_data, data)

  saveRDS(data, RDS_PATH)

  return(RDS_PATH)
}
