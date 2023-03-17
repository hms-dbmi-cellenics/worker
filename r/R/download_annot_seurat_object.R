DownloadAnnotSeuratObject <- function(req, data) {

  cell_sets <- req$body$cellSets

  parsed_cellsets <- parse_cellsets(cell_sets)
  data <- add_clusters(data, parsed_cellsets)

  embedding_data <- req$body$embedding
  data <- assignEmbedding(embedding_data, data)

  fpath <- "/R/r.rds"
  saveRDS(data, fpath)

  return(fpath)
}


assignEmbedding <- function(embedding_data, data, reduction_method = "umap") {
  cells_id <- data@meta.data$cells_id
  embedding <- do.call(rbind, embedding_data)

  # Add 1 to cells_id because it's 0-index and embeddings is not.
  embedding <- embedding[cells_id + 1, ]
  rownames(embedding) <- colnames(data)

  embedding_key = "UMAP_"
  if(reduction_method == "tsne") {
    embedding_key = "tSNE_"
  }

  colnames(embedding) <- c(paste0(embedding_key, "1"), paste0(embedding_key, "2"))
  data[[reduction_method]] <- Seurat::CreateDimReducObject(embeddings = embedding, key = embedding_key)

  return(data)
}

