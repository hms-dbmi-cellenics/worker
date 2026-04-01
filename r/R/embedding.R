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
  use_saved <- req$body$use_saved
  config <- req$body$config
  pca_nPCs <- 30

  set.seed(ULTIMATE_SEED)

  # To run embedding, we need to set the reduction.
  if ("active.reduction" %in% names(data@misc)) {
    active_reduction <- data@misc[["active.reduction"]]
  } else {
    active_reduction <- "pca"
  }

  # The slot numPCs is set in dataIntegration with the selected PCA by the user.
  if ("numPCs" %in% names(data@misc)) {
    pca_npcs <- data@misc[["numPCs"]]
  }

  if (method == "images") {
    img_names <- Seurat::Images(data)
    df_embeddings <- lapply(img_names, get_rotated_tissue_coords, data)
    df_embedding <- do.call(rbind, df_embeddings)

  } else {

    if (!use_saved) {
      data <- getEmbedding(
        config,
        method,
        active_reduction,
        pca_npcs,
        data
      )
    }

    df_embedding <- Seurat::Embeddings(data, reduction = method)
  }

  # Order embedding by cells id in ascending form
  meta <- data@meta.data
  df_embedding <- as.data.frame(df_embedding)
  df_embedding$cells_id <- meta[row.names(df_embedding), "cells_id"]
  df_embedding <- df_embedding[order(df_embedding$cells_id), ]
  df_embedding <- df_embedding |>
    tidyr::complete(cells_id = seq(0, max(meta$cells_id))) |>
    dplyr::select(-"cells_id")

  map2_fun <- function(x, y) {
    if (is.na(x)) {
      NULL
    } else {
      c(x, y)
    }
  }
  res <- purrr::map2(df_embedding[[1]], df_embedding[[2]], map2_fun)
  return(res)
}

# for Visium, data are flipped and rotated before SpatialDimPlot.
# returns the rotated/flipped tissue Coordinates from a Seurat object.
get_rotated_tissue_coords <- function(img_name, scdata) {
  coord_spot <- SeuratObject::GetTissueCoordinates(
    scdata,
    img_name,
    scale = "lowres"
  )[, 2:1] # rotation
  colnames(coord_spot) <- c("x", "y")
  return(coord_spot)
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

  has_sketch <- "sketch" %in% names(data@assays)

  message("Calculating embedding: ",
    "\n- method: ", method,
    "\n- reduction to use: ", reduction_type,
    "\n- number of PCs: ", num_pcs
  )

  if (method == "tsne") {

    # TSNE doesn't support projecting to full data,
    # so we run on full data even if sketch is available
    if (has_sketch) {
      reduction_type <- gsub("[.]sketch$", "", reduction_type)
      warning(
        "TSNE doesn't support sketch projection.",
        "\n- reduction switched to: ", reduction_type, "\n"
      )
    }

    data <- Seurat::RunTSNE(
      data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate
    )

  } else if (method == "umap") {

    if (has_sketch) {
      # use Python umap directly to avoid uwot Docker segfault
      data <- run_sketch_umap(
        data,
        reduction_model = "umap",
        reduction = reduction_type,
        config = config,
        num_pcs = num_pcs
      )
    } else {
      # standard UMAP for non-sketch data
      data <- Seurat::RunUMAP(
        data,
        reduction = reduction_type,
        dims = 1:num_pcs,
        verbose = FALSE,
        min.dist = config$minimumDistance,
        metric = config$distanceMetric,
        umap.method = "umap-learn",
        seed.use = ULTIMATE_SEED
      )
    }
  }
  return(data)
}

# use umap-learn via reticulate to fit sketch and project to full
# avoids uwot Docker segfault (Seurat forces uwot if return.model=TRUE)
run_sketch_umap <- function(
  object, reduction_model, reduction, config, num_pcs
) {

  # get the sketch data (PCA, harmony, etc)
  # use only num_pcs dimensions to match Seurat's RunUMAP dims=1:num_pcs
  sketch_data <- Seurat::Embeddings(object, reduction = reduction)[, 1:num_pcs]

  # import python umap
  umap <- reticulate::import("umap")
  sklearn <- reticulate::import("sklearn.utils")

  # create UMAP model matching Seurat's RunUMAP defaults for umap-learn
  random_state <- sklearn$check_random_state(as.integer(ULTIMATE_SEED))
  umap_model <- umap$UMAP(
    n_neighbors = 15L,
    n_components = 2L,
    metric = config$distanceMetric,
    min_dist = config$minimumDistance,
    random_state = random_state,
    verbose = FALSE
  )

  # fit on sketch and get embedding
  message("Fitting UMAP on downsampled reduction: ", reduction)
  sketch_embedding <- umap_model$fit_transform(as.matrix(sketch_data))

  # set rownames to match cell identifiers from sketch data
  rownames(sketch_embedding) <- rownames(sketch_data)

  # Get full data for projection
  # use only num_pcs dimensions to match Seurat's RunUMAP dims=1:num_pcs
  full_reduction <- gsub("[.]sketch$", "", reduction)
  full_data <- Seurat::Embeddings(object, reduction = full_reduction)
  full_data <- full_data[, 1:num_pcs]

  message("Projecting UMAP to full reduction: ", full_reduction)
  full_embedding <- umap_model$transform(as.matrix(full_data))

  # set rownames to match cell identifiers from full data
  rownames(full_embedding) <- rownames(full_data)

  # store full embedding in main umap reduction
  full_umap_reduction <- Seurat::CreateDimReducObject(
    embeddings = full_embedding,
    key = "UMAP_",
    assay = Seurat::DefaultAssay(object),
    global = TRUE
  )
  object[[reduction_model]] <- full_umap_reduction

  return(object)
}


# assignEmbedding
# Assigns embedding from embedding json to the Seurat object.
# embedding_data is the embedding coordinates.
# data is the seurat object.
#
#' @export
assignEmbedding <- function(embedding_data, data, reduction_method = "umap") {
  cells_id <- data@meta.data$cells_id
  embedding <- do.call(rbind, embedding_data)

  # Add 1 to cells_id because it's 0-index and embeddings is not.
  embedding <- embedding[cells_id + 1, ]
  rownames(embedding) <- colnames(data)

  reduction_keys <- list(
    "umap" = "UMAP_",
    "tsne" = "tSNE_"
  )
  embedding_key <- unname(unlist(reduction_keys[reduction_method]))

  colnames(embedding) <- paste(embedding_key, 1:2, sep = "")

  reduction <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = embedding_key,
    assay = "RNA"
  )
  data[[reduction_method]] <- reduction
  return(data)
}
