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
  pca_npcs <- 30

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
    df_embeddings <- lapply(img_names, function(img_name) {
      scale <- get_image_scale(img_name, data)
      SeuratObject::GetTissueCoordinates(data, img_name, scale = scale)
    })
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
    df_embedding <- as.data.frame(df_embedding[, 1:2])
    colnames(df_embedding) <- c("x", "y")
    df_embedding$cell <- row.names(df_embedding)
  }

  # Order embedding by cells id in ascending form
  meta <- data@meta.data
  df_embedding$cells_id <- meta[df_embedding$cell, "cells_id"]

  df_embedding <- dplyr::arrange(df_embedding, cells_id)

  data.table::setDT(df_embedding)
  result <- vector("list", max(meta$cells_id) + 1L)

  pairs <- df_embedding[,
    .(v = list(c(rbind(x, y)))),
    by = cells_id
  ]

  result[pairs$cells_id + 1L] <- pairs$v
  result
}

get_image_scale <- function(img_name, scdata) {
  # Xenium FOVs have no scale step; Visium HD applies the hires/lowres factor.
  if (is_scaleless_spatial(scdata)) return(NULL)

  dims <- dim(scdata[[img_name]]@image)
  ifelse(any(dims[1:2] <= 600), "lowres", "hires")
}

# Spatial technologies whose coordinates have no pixel->coord scale step (the
# accessors return micron coords directly). Dispatch on the technology persisted
# onto the object during create-seurat (data@misc$technology), not on object
# introspection: the worker loads the saved Seurat object, not the pipeline
# config, so the persisted technology is the authoritative signal.
SCALELESS_SPATIAL_TECHNOLOGIES <- c("xenium")

is_scaleless_spatial <- function(scdata) {
  technology <- scdata@misc$technology
  !is.null(technology) && technology %in% SCALELESS_SPATIAL_TECHNOLOGIES
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

    data <- run_umap(
      data,
      reduction_model = "umap",
      reduction = reduction_type,
      config = config,
      num_pcs = num_pcs,
      has_sketch = has_sketch
    )
  }

  return(data)
}

# Number of neighbours used to build the UMAP graph. uwot's default is 15;
# 30 is Seurat's, which is what the platform has always used.
UMAP_N_NEIGHBORS <- 30L

# HNSW index parameters, matching the defaults uwot::umap2 uses internally
# (uwot:::hnsw_nn) so that computing the neighbours here is only a change of
# who runs the build, not of what gets built.
HNSW_M <- 16
HNSW_EF_CONSTRUCTION <- 200
HNSW_EF <- 10

# uwot builds its HNSW index with as many threads as it is given, and a
# multithreaded build is not reproducible: the same input gives a different
# neighbour graph on each run. Building single threaded is the one lever that
# makes it deterministic (~10s at 47k cells x 28 dims, against ~2s
# multithreaded), and it is a build-only cost: search stays multithreaded.
build_umap_nn_index <- function(x, metric) {
  RcppHNSW::hnsw_build(
    x,
    distance = hnsw_distance(metric),
    M = HNSW_M,
    ef = HNSW_EF_CONSTRUCTION,
    n_threads = 1,
    verbose = FALSE
  )
}

# Nearest neighbours of x in an index, in the list form uwot accepts for
# nn_method. Used both to fit (x is the fitted data) and to project a sketch
# model onto the full data (x is the full data, index built on the sketch).
search_umap_nn_index <- function(x, index, metric) {
  res <- RcppHNSW::hnsw_search(
    X = x,
    k = UMAP_N_NEIGHBORS,
    ann = index,
    ef = HNSW_EF,
    n_threads = parallel::detectCores(),
    verbose = FALSE
  )

  # uwot indexes euclidean data with the squared L2 metric, so undo it here too
  if (metric == "euclidean") res$dist <- sqrt(res$dist)

  list(idx = res$idx, dist = res$dist)
}

# RcppHNSW has no euclidean class, uwot uses l2 and square roots the distances
hnsw_distance <- function(metric) {
  if (metric == "euclidean") "l2" else metric
}

# uwot implementation of UMAP, and the only one that can fit on a sketch and
# project onto the full dataset (Seurat forces uwot when return.model = TRUE).
run_umap <- function(
  object, reduction_model, reduction, config, num_pcs, has_sketch = FALSE
) {

  red_data <- as.matrix(
    Seurat::Embeddings(object, reduction = reduction)[, 1:num_pcs]
  )
  metric <- config$distanceMetric

  tstart_umap <- Sys.time()
  message("Fitting UMAP on ", reduction, " data via uwot...")

  # NOTE: uwot::umap2 with RcppHNSW installed used as faster and
  # avoids segfaults seen from RcppAnnoy with n_threads != 0
  #
  # Reproducibility: the seed, batch = TRUE, fast_sgd = FALSE and
  # rng_type = "deterministic" (uwot >= 0.2.3) cover the optimisation phase, but
  # they are not enough on their own: uwot builds its HNSW index with n_threads
  # and a multithreaded build gives a different neighbour graph - and so a
  # different layout - on every run. uwot's own n_build_threads is not in CRAN
  # 0.2.4, so we build the index single threaded ourselves and hand uwot the
  # neighbours. Index search is deterministic, so it keeps every thread.
  nn_index <- build_umap_nn_index(red_data, metric)

  set.seed(ULTIMATE_SEED)
  umap_res <- uwot::umap2(
    X = red_data,
    nn_method = search_umap_nn_index(red_data, nn_index, metric),
    n_neighbors = UMAP_N_NEIGHBORS,
    min_dist = config$minimumDistance,
    metric = metric,
    ret_model = has_sketch,
    n_threads = parallel::detectCores(),
    n_sgd_threads = "auto",
    batch = TRUE,
    fast_sgd = FALSE,
    rng_type = "deterministic",
    seed = as.integer(ULTIMATE_SEED)
  )

  message(
    "UMAP fit time: ",
    round(difftime(Sys.time(), tstart_umap, units = "secs"), 2), " seconds"
  )

  if (!has_sketch) {
    full_embedding <- umap_res

  } else {
    full_data <- as.matrix(
      Seurat::Embeddings(
        object,
        reduction = gsub("[.]sketch$", "", reduction)
      )[, 1:num_pcs]
    )

    message("Projecting full dataset...")
    tstart_project <- Sys.time()

    # the model carries no index (it was fit from neighbours we computed), so
    # the full data is searched against the same single threaded index here.
    # rng_type comes from the model; batch/seed still have to be passed.
    full_embedding <- uwot::umap_transform(
      nn_method = search_umap_nn_index(full_data, nn_index, metric),
      model = umap_res,
      n_threads = parallel::detectCores(),
      n_sgd_threads = "auto",
      batch = TRUE,
      seed = as.integer(ULTIMATE_SEED)
    )

    # projecting from neighbours rather than from X means uwot has no cell
    # names to carry over, and CreateDimReducObject requires them
    rownames(full_embedding) <- rownames(full_data)

    message(
      "UMAP projection time: ",
      round(difftime(Sys.time(), tstart_project, units = "secs"), 2), " seconds"
    )
  }

  colnames(full_embedding) <- paste0("UMAP_", 1:2)

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

  # Add 1 to cells_id because it's 0-index and embeddings is not.
  embedding_data <- embedding_data[cells_id + 1]
  embedding <- do.call(rbind, embedding_data)

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
