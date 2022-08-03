# runEmbedding
# Returns a list of x,y coordinates for each cell.
# req is the request.
#
# Req$body has:
# type = type of embedding, supported umap, pca and tsne.
# config = config list.
#
# Config has=
# UMAP:
# minimumDistance = float
# distanceMetric = string (euclidean, cosine, etc)
#
# tsne:
# perplexity
# lerarningRate
#
#' @export
runEmbedding <- function(req, data) {
  method <- req$body$type
  config <- req$body$config
  pca_nPCs <- 30

  set.seed(ULTIMATE_SEED)

  # To run embedding, we need to set the reduction.
  if ("active.reduction" %in% names(data@misc)) {
    active.reduction <- data@misc[["active.reduction"]]
  } else {
    active.reduction <- "pca"
  }

  # The slot numPCs is set in dataIntegration with the selected PCA by the user.
  if ("numPCs" %in% names(data@misc)) {
    pca_nPCs <- data@misc[["numPCs"]]
  }

  message("Active reduction --> ", active.reduction)
  message("Active numPCs --> ", pca_nPCs)
  message("Number of cells/sample:")
  table(data$samples)

  data <- getEmbedding(config, method, active.reduction, pca_nPCs, data)
  df_embedding <- Seurat::Embeddings(data, reduction = method)

  # Order embedding by cells id in ascending form
  df_embedding <- as.data.frame(df_embedding)
  df_embedding$cells_id <- data@meta.data$cells_id
  df_embedding <- df_embedding[order(df_embedding$cells_id), ]
  df_embedding <- df_embedding %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    dplyr::select(-cells_id)

  map2_fun <- function(x, y) {
    if (is.na(x)) {
      return(NULL)
    } else {
      return(c(x, y))
    }
  }
  res <- purrr::map2(df_embedding[[1]], df_embedding[[2]], map2_fun)
  return(res)
}


# getEmbedding
# Return embedding calculated for the seurat object.
# config is the embedding work request config.
# data is the seurat object.
# method is the embedding method, e.g. UMAP.
# reduction_type is the type of reduction that is used, e.g. PCA.
# num_pcs is the number of principal components.
#
#' @export
getEmbedding <- function(config, method, reduction_type, num_pcs, data) {
  if (method == "tsne") {
    data <- Seurat::RunTSNE(data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate
    )
  } else if (method == "umap") {
    data <- Seurat::RunUMAP(data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      verbose = FALSE,
      min.dist = config$minimumDistance,
      metric = config$distanceMetric,
      umap.method = "umap-learn"
    )
  }

  return(data)
}


# assignEmbedding
# Assigns embedding from embedding json to the Seurat object.
# embedding_data is the embedding coordinates.
# data is the seurat object.
#
#' @export
assignEmbedding <- function(embedding_data, data) {
  # Filter for nulls
  embedding <- embedding_data[!vapply(embedding_data, is.null, logical(1))]
  embedding <- do.call(rbind, embedding)
  rownames(embedding) <- colnames(data)

  data[["umap"]] <- Seurat::CreateDimReducObject(embeddings = embedding, key = "UMAP_")

  return(data)
}
