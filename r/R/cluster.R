# IMPORTANT: functions in this file are duplicated in the pipeline.
# If you update, change both.

#' Get Clusters
#'
#' @param data SeuratObject
#' @param req list with items:
#'  \code{req$body$type} either "louvain" or "leiden"
#'  \code{req$body$config$resolution} integer, range: 0 - 2
#'
#' @return data.frame with columns "cluster" and "cell_ids"
#' @export
#'
#' @examples
runClusters <- function(req, data) {
  type <- req$body$type
  config <- req$body$config
  resolution <- config$resolution

  data <- getClusters(type, resolution, data)
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
  # In the meta data slot the clustering is stored with the resolution
  # used to calculate it
  # RNA_snn_res.#resolution
  meta <- data@meta.data
  df <- data.frame(
    cluster = meta[, res_col],
    cell_ids = meta$cells_id,
    row.names = rownames(meta)
  )

  message("formatting cellsets")
  # for spatial data, use Spaco to assign spatially-coherent cluster colors.
  # non-spatial data (or Spaco failure) keeps the default color pool.
  color_pool <- data@misc$color_pool
  if (length(Seurat::Images(data)) > 0) {
    spaco_colors <- get_spaco_color_map(data, df)
    if (!is.null(spaco_colors)) color_pool <- spaco_colors
  }
  formatted_cell_sets <- format_cell_sets_object(df, type, color_pool)

  message("updating through api")
  updateCellSetsThroughApi(
    formatted_cell_sets,
    req$body$apiUrl,
    req$body$experimentId,
    "louvain",
    req$body$authJwt,
    append = FALSE
  )

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

  # use the reduction from data integration for nearest neighbors graph
  if ("active.reduction" %in% names(data@misc)) {
    active_reduction <- data@misc[["active.reduction"]]
  } else {
    active_reduction <- "pca"
  }

  # run clustering on sketch assay if it exists
  has.sketch <- "sketch" %in% names(data@assays)
  if (has.sketch) Seurat::DefaultAssay(data) <- "sketch"

  res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
  algorithm <- list("louvain" = 1, "leiden" = 4)[[type]]

  if (type == "leiden") {
    # emulate FindClusters: overwrites seurat_clusters slot and meta.data column
    snn_graph <- getSNNiGraph(data, active_reduction)
    clus_res <- igraph::cluster_leiden(
      snn_graph,
      "modularity",
      resolution_parameter = resolution
    )

    clusters <- clus_res$membership
    names(clusters) <- clus_res$names
    clusters <- clusters[colnames(data)]
    data$seurat_clusters <- data@meta.data[, res_col] <- factor(clusters - 1)

  } else {

    graph_name <- paste0(Seurat::DefaultAssay(data), "_snn")
    if (!graph_name %in% names(data)) {
      # number of dimensions used must be lte to available dimensions
      dims <- 1:min(10, length(data@reductions[[active_reduction]]))
      data <-
        Seurat::FindNeighbors(
          data,
          annoy.metric = "cosine",
          reduction = active_reduction,
          dims = dims,
          verbose = FALSE,
        )
    }
    data <- Seurat::FindClusters(
      data,
      resolution = resolution,
      verbose = FALSE,
      algorithm = algorithm,
      random.seed = ULTIMATE_SEED
    )
  }

  if (has.sketch) {
    # prepare to project clusters from sketch to full data
    data$seurat_clusters_sketch <- data$seurat_clusters
    data$seurat_clusters <- NULL

    # IMPORTANT: need to use same dims as previous call to ProjectData
    # otherwise will hang indefinitely
    full_reduction <- gsub("[.]sketch$", "", active_reduction)
    full_npcs <- length(data@reductions[[full_reduction]])

    # extend results to full dataset
    # NOTE: uses stored calculations from previous call (during integration)
    data <- Seurat::ProjectData(
      object = data,
      assay = "RNA",
      full.reduction = full_reduction,
      sketched.assay = "sketch",
      sketched.reduction = active_reduction,
      dims = 1:full_npcs,
      verbose = TRUE,
      refdata = list(seurat_clusters = "seurat_clusters_sketch")
    )

    # transfer res_col from sketch to RNA assay
    # this is what is uploaded to cellsets json for worker
    # also used by worker to check for existing clustering if change resolution
    rna_res_col <- gsub("sketch", "RNA", res_col)
    data@meta.data[, rna_res_col] <- data$seurat_clusters

    # set default assay back to RNA so downstream code (runClusters) reads the
    # full projected clusters (all cells), not the sketch column (sketched only)
    Seurat::DefaultAssay(data) <- "RNA"
    message("Projected clusters from sketch to full data")
  }

  return(data)
}


#' Get and Convert SNN Graph object into igraph object
#'
#' This is used to facilitate leiden clustering.
#'
#' @param data \code{Seurat} object
#'
#' @return boolean indicating if SNN Graph object exists
#'
getSNNiGraph <- function(data, active_reduction) {
  # check to see if we already have Seurat SNN Graph object
  snn_name <- paste0(data@active.assay, "_snn")

  # if doesn't exist, run SNN
  if (!snn_name %in% names(data)) {
    dims <- 1:min(10, length(data@reductions[[active_reduction]]))
    data <-
      Seurat::FindNeighbors(data,
                            reduction = active_reduction,
                            dims = dims,
                            verbose = FALSE)
  }

  # convert Seurat Graph object to igraph
  # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
  adj_matrix <-
    Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE
  )
  return(graph)
}

#' Formats cell sets object for patching through the API
#'
#' This function is only used to format clustering cellsets. Converting from
#' data.frame to list and adding slots necessary for the cellsets file.
#'
#' @param cell_sets data.frame with two columns: cluster and cell_ids
#' @param clustering_method string Either louvain or leiden.
#' @param color_pool character vector of colors in hex, consumed one per cluster
#'  in the order clusters are sorted. For spatial data this is the Spaco-derived
#'  palette (see \code{\link{get_spaco_color_map}}), already ordered to match.
#'
#' @return list
#' @export
#'
#' @examples
format_cell_sets_object <-
  function(cell_sets, clustering_method, color_pool) {
    name <- paste0(clustering_method, " clusters")
    cell_sets_object <-
      list(
        key = "louvain",
        name = name,
        rootNode = TRUE,
        type = "cellSets",
        children = list()
      )
    for (i in sort(unique(cell_sets$cluster))) {
      cells <- cell_sets[cell_sets$cluster == i, "cell_ids"]
      new_set <- list(
        key = paste0(clustering_method, "-", i),
        name = paste0("Cluster ", i),
        rootNode = FALSE,
        type = "cellSets",
        color = color_pool[1],
        cellIds = unname(cells)
      )
      color_pool <- color_pool[-1]
      cell_sets_object$children <-
        append(cell_sets_object$children, list(new_set))
    }
    return(cell_sets_object)
  }

#' Compute spatially-coherent cluster colors (Spaco algorithm, in R)
#'
#' For spatial datasets (e.g. Visium HD), assigns cluster colors so that
#' spatially interlaced clusters get perceptually distinct colors. This is a
#' native R reimplementation of the Spaco algorithm
#' (\url{https://github.com/BrainStOrmics/Spaco}); each image in the Seurat
#' object is treated as a slice. We reimplement it rather than calling the
#' Python package to avoid its slow numba start-up.
#'
#' The pipeline matches Spaco: per slice, build a cluster interlacement distance
#' graph (\code{\link{spatial_distance_r}}), sum the graphs across slices, then
#' assign \code{SPACO_COLOR_POOL} so the most interlaced clusters get the most
#' distinct colors (\code{\link{assign_spaco_palette}}).
#'
#' When the data is sketched, colors are computed from the sketch cells only
#' (far fewer cells, so the neighbourhood graph is much cheaper); cluster colors
#' apply to the whole dataset regardless.
#'
#' @param scdata Seurat object with one or more images
#' @param cell_sets data.frame with columns \code{cluster} and \code{cell_ids},
#'  with cell barcodes as row names.
#'
#' @return character vector of hex colors, ordered to match the sorted cluster
#'  order used by \code{\link{format_cell_sets_object}}, or \code{NULL} if no
#'  slice has cluster-labelled cells or anything fails.
#' @export
#'

# Colors Spaco assigns to clusters (mirrors the pipeline's SPACO_COLOR_POOL),
# most distinct first: the first 20 of the
# default color pool, then colors each chosen as the one most different
# (CIEDE2000) from all before it (as in Glasbey et al. 2007), up to 100 colors.
# Candidates are an sRGB grid of 18 levels per channel excluding dark (Lab L <= 30)
# and grey (chroma <= 15) colors. Colors within CIEDE2000 25 of the UI's light
# grey background (#F0F2F5), where they are hard to see, are excluded throughout.
# The default pool alone has near-identical colors after its first ~20, which
# interlaced clusters of large datasets would get.
SPACO_COLOR_POOL <- c(
  "#77aadd", "#ee8866", "#eedd88", "#ffaabb", "#44bb99", "#bbcc33",
  "#999900", "#bc9cc9", "#295a8b", "#b14d2c", "#c2aa36", "#d4788a",
  "#3e8ab1", "#158364", "#737a38", "#4e4e00", "#835794", "#edf761",
  "#4d8cc4", "#96003c", "#ff00ff", "#00ff69", "#4b00e1", "#00a50f",
  "#1ee1f0", "#c37800", "#8769ff", "#2d5a5a", "#693c0f", "#006900",
  "#e10087", "#967869", "#a5b487", "#ffb469", "#ff003c", "#78ffd2",
  "#965a69", "#870078", "#b400e1", "#78a5a5", "#787896", "#96000f",
  "#a59678", "#4b3c87", "#ff96ff", "#876900", "#c369b4", "#693c3c",
  "#005a3c", "#ff3c00", "#a5e196", "#b496ff", "#78694b", "#005aff",
  "#699669", "#ff7800", "#69d200", "#ffc300", "#3c8787", "#00b4d2",
  "#d2004b", "#965a1e", "#4b5a3c", "#5a4b69", "#2dc369", "#785a4b",
  "#a5873c", "#c39678", "#ff69b4", "#9678c3", "#a5ff00", "#6978c3",
  "#783c5a", "#c3695a", "#2d873c", "#784bc3", "#ff7878", "#c39696",
  "#ff3c78", "#006987", "#87a54b", "#5a4b2d", "#8796b4", "#967896",
  "#00b4b4", "#5a785a", "#c3692d", "#d29600", "#b43ca5", "#d2001e",
  "#003cc3", "#b45ae1", "#00f0e1", "#3cc3ff", "#8700b4", "#e1ffa5",
  "#963c3c", "#872d00", "#a5784b", "#5a8700"
)


# Spaco neighbourhood radius, in microns (the Spaco paper's default). Applied
# directly to Xenium (micron coords); converted to image pixels for Visium HD.
SPACO_RADIUS_MICRONS <- 50

#' Neighbourhood radius for \code{spatial_distance_r}, in the coordinate units
#' returned by \code{GetTissueCoordinates(scdata, image_name, scale = scale)}.
#'
#' Mirrors the pipeline's \code{get_spaco_radius}. The Spaco radius is defined in
#' microns. Xenium coords are already microns (\code{scale} is NULL) so it's used
#' as-is; Visium HD coords are image pixels at the hires/lowres scale, so convert
#' microns -> fullres pixels (via \code{microns_per_pixel}, persisted onto
#' \code{@misc} by the pipeline) -> scaled pixels (via the image scale factor).
#' Returns \code{NULL} when the factor is unavailable so the caller skips the
#' slice rather than running the neighbourhood graph at the wrong scale.
get_spaco_radius <- function(scdata, image_name, scale) {
  if (is.null(scale)) {
    return(SPACO_RADIUS_MICRONS)
  }

  microns_per_pixel <- scdata@misc$microns_per_pixel
  scale_factor <- switch(scale,
    hires = scdata[[image_name]]@scale.factors$hires,
    lowres = scdata[[image_name]]@scale.factors$lowres,
    NULL
  )
  if (is.null(microns_per_pixel) || is.null(scale_factor)) {
    return(NULL)
  }

  SPACO_RADIUS_MICRONS / microns_per_pixel * scale_factor
}

get_spaco_color_map <- function(scdata, cell_sets) {
  tstart <- Sys.time()
  tryCatch(
    {
      img_names <- Seurat::Images(scdata)
      if (length(img_names) == 0) {
        return(NULL)
      }

      # if sketched, compute colors from the sketch cells only (far fewer);
      # the resulting cluster colors apply to the full dataset
      cells_use <- NULL
      if ("sketch" %in% names(scdata@assays)) {
        cells_use <- colnames(scdata[["sketch"]])
      }

      clusters <- as.character(sort(unique(cell_sets$cluster)))

      # per-slice cluster interlacement distance graphs
      per_slice <- lapply(img_names, function(img) {
        scale <- get_image_scale(img, scdata)
        coords <- SeuratObject::GetTissueCoordinates(scdata, img, scale = scale)
        labels <- as.character(cell_sets[coords$cell, "cluster"])
        keep <- !is.na(labels)
        if (!is.null(cells_use)) keep <- keep & coords$cell %in% cells_use
        if (sum(keep) < 4) {
          return(NULL)
        }
        # 50um neighbourhood in this slice's coordinate units; skip the slice if
        # the micron->pixel factor is missing rather than colour at a wrong scale
        radius <- get_spaco_radius(scdata, img, scale)
        if (is.null(radius)) {
          return(NULL)
        }
        spatial_distance_r(
          as.matrix(coords[keep, c("x", "y")]), labels[keep], radius = radius
        )
      })
      per_slice <- Filter(Negate(is.null), per_slice)
      if (length(per_slice) == 0) {
        return(NULL)
      }

      # sum interlacement graphs across slices, aligned to the full cluster set
      cluster_distance <- merge_cluster_distances(per_slice, clusters)

      # most interlaced clusters get the most distinct colors
      color_map <- assign_spaco_palette(cluster_distance, SPACO_COLOR_POOL)

      message(
        "Spaco(R): coloring time ",
        round(difftime(Sys.time(), tstart, units = "secs"), 2),
        " seconds for ", length(clusters), " clusters"
      )

      palette <- unname(color_map[clusters])
      # every cluster should get a colour; if any is missing (e.g. a cluster
      # absent from the interlacement graph) fall back rather than hand an NA
      # colour to format_cell_sets_object (mirrors the pipeline's spaco_palette)
      if (anyNA(palette)) {
        return(NULL)
      }
      palette
    },
    error = function(e) {
      message(
        "Spaco(R) coloring failed (", conditionMessage(e),
        "), using default colors"
      )
      NULL
    }
  )
}

#' Cluster spatial interlacement distance graph (Spaco \code{spatial_distance})
#'
#' Native R port of \code{spaco.distance.spatial_distance}. For each cell, finds
#' its \code{n_neighbors} nearest neighbours within \code{radius}; cells with
#' fewer than \code{n_cells} same-cluster neighbours are "banished" (contribute
#' nothing). The interlacement score between two clusters accumulates
#' inverse-distance weights over neighbouring cells of different clusters, then
#' is symmetrised (max) with a zero diagonal.
#'
#' @param coords numeric matrix of spatial coordinates (n cells x 2)
#' @param labels character vector of cluster labels (length n)
#' @param radius,n_neighbors,n_cells Spaco neighbourhood parameters. \code{radius}
#'  is in the same units as \code{coords}; production callers always pass a value
#'  derived from \code{\link{get_spaco_radius}} (50um converted to the coords'
#'  scale). The default is only a fallback for direct/unit-test calls.
#'
#' @return symmetric cluster x cluster matrix with cluster labels as dimnames
#'
spatial_distance_r <- function(
  coords, labels, radius = 90, n_neighbors = 16L, n_cells = 3L
) {
  clusters <- sort(unique(labels))
  k_clusters <- length(clusters)
  n <- nrow(coords)

  knn <- RANN::nn2(coords, coords, k = min(n_neighbors, n))
  idx <- knn$nn.idx     # n x k, self in column 1
  dst <- knn$nn.dists   # n x k, self distance 0 in column 1

  label_idx <- match(labels, clusters)
  neighbor_label <- matrix(label_idx[idx], nrow = n)
  within <- dst <= radius

  # banish cells with too few same-cluster neighbours within radius
  same_count <- rowSums((neighbor_label == label_idx) & within)
  banished <- same_count < n_cells
  # neighbourhood size (including self) used to normalise scores
  size_n <- rowSums(within)

  # accumulate scores over non-self, within-radius neighbours of kept cells
  mask <- within
  mask[, 1] <- FALSE
  mask[banished, ] <- FALSE
  pos <- which(mask, arr.ind = TRUE)

  weights <- (1 / dst[pos]) / size_n[pos[, 1]]
  cluster_i <- label_idx[pos[, 1]]
  cluster_j <- neighbor_label[pos]
  agg <- tapply(
    weights,
    list(factor(cluster_i, seq_len(k_clusters)),
         factor(cluster_j, seq_len(k_clusters))),
    sum
  )
  agg[is.na(agg)] <- 0
  score <- matrix(as.numeric(agg), k_clusters, k_clusters)

  # keep max of the two directions, zero the diagonal
  score <- pmax(score, t(score))
  diag(score) <- 0
  dimnames(score) <- list(clusters, clusters)
  score
}

#' Sum per-slice cluster interlacement graphs, aligned to a cluster set
#'
#' @param per_slice list of cluster x cluster matrices (one per slice)
#' @param clusters character vector of all clusters to align to
#'
#' @return cluster x cluster matrix summed across slices (missing entries 0)
#'
merge_cluster_distances <- function(per_slice, clusters) {
  merged <- matrix(
    0, length(clusters), length(clusters), dimnames = list(clusters, clusters)
  )
  for (m in per_slice) {
    rn <- rownames(m)
    merged[rn, rn] <- merged[rn, rn] + m
  }
  merged
}

#' Assign palette colors to clusters by interlacement (Spaco palette mode)
#'
#' Mirrors the pipeline's \code{assign_spaco_palette}. Spaco's mode for a given
#' palette: choose which cluster gets which color so
#' that the most interlaced clusters get the most distinct colors, maximising
#' the sum over cluster pairs of interlacement x CIEDE2000 color difference.
#' Clusters are assigned greedily, most interlaced first, then pairs of clusters
#' swap colors while that improves the sum.
#'
#' Used instead of Spaco's auto-generated palette (a UMAP embedding of the graph
#' into CIE-Lab), which crowds many clusters into similar hues and changes from
#' run to run.
#'
#' @param cluster_distance symmetric cluster x cluster interlacement matrix
#' @param palette hex colors, most distinct first. The first one per cluster are
#'  used, recycled if there are more clusters than colors.
#' @param max_passes maximum passes over the clusters when swapping. Almost all
#'  of the improvement comes in the first few, and each costs ~k^3 operations
#'
#' @return named character vector mapping cluster to hex color
#'
assign_spaco_palette <- function(cluster_distance, palette, max_passes = 3) {
  clusters <- rownames(cluster_distance)
  k <- length(clusters)
  colors <- rep_len(palette, k)

  lab <- farver::convert_colour(t(grDevices::col2rgb(colors)), "rgb", "lab")
  color_diff <- farver::compare_colour(lab, lab, "lab", method = "cie2000")

  # assigned[i] is the index in colors of cluster i's color. Ties go to the
  # earlier color, so unconnected clusters take the palette in order.
  assigned <- rep(NA_integer_, k)
  for (i in order(-rowSums(cluster_distance))) {
    done <- which(!is.na(assigned))
    free <- setdiff(seq_len(k), assigned[done])
    gain <- vapply(free, function(color) {
      sum(cluster_distance[i, done] * color_diff[color, assigned[done]])
    }, numeric(1))
    assigned[i] <- free[which.max(gain)]
  }

  # swap cluster i's color with the cluster j where that improves the sum most.
  # Swapping only changes their pairs with the other clusters m:
  # gain_j = sum_m (d[i, m] - d[j, m]) * (c[j, m] - c[i, m]), m != i, j
  for (pass in seq_len(max_passes)) {
    improved <- FALSE
    for (i in seq_len(k)) {
      pair_diff <- color_diff[assigned, assigned]
      gain <- rowSums(
        sweep(-cluster_distance, 2, cluster_distance[i, ], `+`) *
          sweep(pair_diff, 2, pair_diff[i, ])
      ) + 2 * cluster_distance[i, ] * pair_diff[i, ] # remove the m = i, j terms
      gain[i] <- 0

      j <- which.max(gain)
      if (gain[j] > 1e-9) {
        assigned[c(i, j)] <- assigned[c(j, i)]
        improved <- TRUE
      }
    }
    if (!improved) break
  }

  stats::setNames(colors[assigned], clusters)
}
