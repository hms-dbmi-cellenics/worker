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
  print(table(data$samples))

  if (method == 'images') {
    img_names <- Seurat::Images(data)
    df_embeddings <- lapply(img_names, get_rotated_tissue_coords, data)
    df_embedding <- do.call(rbind, df_embeddings)

  } else {
    if (!use_saved)
      data <- getEmbedding(config, method, active.reduction, pca_nPCs, data)

    df_embedding <- Seurat::Embeddings(data, reduction = method)
  }

  # Order embedding by cells id in ascending form
  meta <- data@meta.data
  df_embedding <- as.data.frame(df_embedding)
  df_embedding$cells_id <- meta[row.names(df_embedding), "cells_id"]
  df_embedding <- df_embedding[order(df_embedding$cells_id), ]
  df_embedding <- df_embedding |>
    tidyr::complete(cells_id = seq(0, max(meta$cells_id))) |>
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

# at least in the case of Visium, data are flipped and rotated before SpatialDimPlot.
# this function  return the rotated/flipped tissue Coordinates from a Seurat object.
get_rotated_tissue_coords <- function(img_name, scdata) {
  coord_spot <- SeuratObject::GetTissueCoordinates(scdata, img_name, scale = "lowres")[,2:1] # rotation
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

  if (method == "tsne") {

    data <- Seurat::RunTSNE(
      data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate
    )

  } else if (method == "umap") {

    has_sketch <- "sketch" %in% names(data@assays)
    if (has_sketch) {
      # Use Python umap directly to avoid uwot Docker segfault
      data <- runSketchUMAP(
        data,
        reduction.model = "umap",
        reduction = reduction_type,
        config = config,
        num_pcs = num_pcs
      )
    } else {
      # Standard UMAP for non-sketch data
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
runSketchUMAP <- function(object, reduction.model, reduction, config, num_pcs) {

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
    verbose = TRUE
  )

  # fit on sketch and get embedding
  message("Fitting UMAP on sketch assay")
  sketch_embedding <- umap_model$fit_transform(as.matrix(sketch_data))

  # set rownames to match cell identifiers from sketch data
  rownames(sketch_embedding) <- rownames(sketch_data)

  # store sketch embedding in umap.sketch reduction
  reduction.sketch <- paste0(reduction.model, ".sketch")
  sketch_umap_reduction <- Seurat::CreateDimReducObject(
    embeddings = sketch_embedding,
    key = "UMAP.sketch_",
    assay = Seurat::DefaultAssay(object[[reduction]]),
    global = TRUE
  )
  object[[reduction.sketch]] <- sketch_umap_reduction

  # Get full data for projection
  # use only num_pcs dimensions to match Seurat's RunUMAP dims=1:num_pcs
  full_reduction <- gsub("[.]sketch$", "", reduction)
  full_data <- Seurat::Embeddings(object, reduction = full_reduction)
  full_data <- full_data[, 1:num_pcs]

  message("Projecting UMAP embedding from sketch to full dataset")
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
  object[[reduction.model]] <- full_umap_reduction

  message("Done projecting UMAP embedding from sketch to full dataset")
  return(object)
}

# Fallback: needed because ProjectData doesn't support umap-learn
# uwot segfaults on Docker
ProjectSketchedUMAP <- function(object, reduction.model, reduction) {
  library(Seurat)

  full_sketch.nn <- Tool(object = object, slot = "TransferSketchLabels")$full_sketch.nn
  full_sketch.weight <- Tool(object = object, slot = "TransferSketchLabels")$full_sketch.weight

  umap.model <- Misc(object = object[[reduction.model]], slot = "model")

  if (ncol(full_sketch.nn) > umap.model$n_neighbors) {
    full_sketch.nn@nn.idx <-
      full_sketch.nn@nn.idx[, 1:umap.model$n_neighbors]

    full_sketch.nn@nn.dist <-
      full_sketch.nn@nn.dist[, 1:umap.model$n_neighbors]
  }

  message("Running UMAP projection from sketch to full dataset")

  proj.umap <- RunUMAP(
    object = full_sketch.nn,
    reduction.model = object[[reduction.model]],
    verbose = TRUE,
    assay = slot(object = object[[reduction]], name = "assay.used"),
    umap.method = "umap-learn"
  )

  # move sketched umap to umap.sketch reduction
  reduction.sketch <- paste0(reduction.model, ".sketch")
  umap.sketch <- object[[reduction.model]]
  Key(umap.sketch) <- Key(reduction.sketch)
  object[[reduction.sketch]] <- umap.sketch

  # set projected umap to umap reduction
  object[[reduction.model]] <- proj.umap
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
  data[[reduction_method]] <- Seurat::CreateDimReducObject(embeddings = embedding, key = embedding_key, assay = "RNA")

  return(data)
}
