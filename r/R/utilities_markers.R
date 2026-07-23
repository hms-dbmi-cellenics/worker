#' getTopMarkerGenes
#'
#' Uses presto::wilcoxauc to find the marker genes that distinguish the
#' cellsets. It then filters the list of genes up to nfeatures using reasonable
#' defaults.
#'
#' @param nfeatures int number of marker genes to get
#' @param data SeuratObject
#' @param cell_sets_ids list of cell sets to split for marker gene selection
#' @param auc_min min area under the wilcoxon test's ROC for a gene to be
#'  considered a marker
#' @param pct_in_min min percentage of cells in cellset that have to express a
#'  gene for it to be considered a marker
#' @param pct_out_max max percentage of cells outside cellset that can express a
#'  gene for it to be considered a marker
#'
#' @return data.frame of top marker genes
#' @export
#'
getTopMarkerGenes <- function(
  nfeatures,
  data,
  cell_sets_ids,
  auc_min = 0.3,
  pct_in_min = 20,
  pct_out_max = 70
) {
  all_markers <- computeMarkerStats(data, cell_sets_ids)$markers

  # may not return nfeatures markers per cluster if values are too stringent
  filtered_markers <- all_markers |>
    dplyr::filter(logFC > 0 &
        auc >= auc_min &
        pct_in >= pct_in_min &
        pct_out <= pct_out_max
    ) |>
    dplyr::arrange(pval) |>
    dplyr::distinct(feature, .keep_all = TRUE) |>
    dplyr::arrange(-logFC)

  # Get top nfeatures per group without additional sorting
  top_markers <- filtered_markers |>
    dplyr::group_by(group) |>
    dplyr::slice_head(n = nfeatures) |>
    dplyr::ungroup() |>
    dplyr::arrange(group)

  return(top_markers)
}

#' computeMarkerStats
#'
#' Factored out of [getTopMarkerGenes()]: downsamples the cells (for speed) and
#' runs presto::wilcoxauc to compute one-vs-rest statistics for every gene in
#' each group. Returns the raw wilcoxauc table together with the subsetted count
#' matrix and per-cell group labels, so callers can compute additional
#' statistics (e.g. Seurat-style fold changes) on the exact same cells.
#'
#' Groups are numbered 1..length(cell_sets_ids) in the order of `cell_sets_ids`.
#'
#' @param data SeuratObject
#' @param cell_sets_ids list of cell sets (vectors of cell ids) to split for
#'  marker gene selection
#' @param downsample_n int max number of cells sampled per group
#' @param seed int RNG seed for the downsampling
#'
#' @return list with:
#'  \code{markers} data.frame from presto::wilcoxauc (group as numeric),
#'  \code{mat} dgCMatrix of the subsetted `data` layer (genes x sampled cells),
#'  \code{groups} numeric vector of group labels for the sampled cells
#' @export
#'
computeMarkerStats <- function(
  data,
  cell_sets_ids,
  downsample_n = 1000,
  seed = 0
) {
  object_ids <- data$cells_id

  # build cell-to-group mapping using data.table for speed
  cell_group_map <- data.table::rbindlist(
    lapply(seq_along(cell_sets_ids), function(i) {
      data.table::data.table(
        cell_id = unlist(cell_sets_ids[[i]]),
        group = i
      )
    })
  )

  # join to assign marker groups to all cells
  dt <- data.table::data.table(cell_id = object_ids)
  dt <- cell_group_map[dt, on = "cell_id"]
  marker_groups <- dt$group

  # for speed: take at most downsample_n cells per cluster
  set.seed(seed)

  # Create temp data frame for sampling
  keep_cell_ids <- data@meta.data |>
    dplyr::mutate(marker_groups = marker_groups) |>
    dplyr::group_by(marker_groups) |>
    dplyr::slice_sample(n = downsample_n, replace = FALSE) |>
    dplyr::ungroup() |>
    dplyr::pull(cells_id)

  # Extract and subset matrix
  mat <- data[["RNA"]]$data
  keep_indices <- match(keep_cell_ids, object_ids)
  marker_groups_subset <- marker_groups[keep_indices]
  mat_subset <- mat[, keep_indices]
  mat_subset <- as(mat_subset, "dgCMatrix")

  all_markers <- presto::wilcoxauc(
    mat_subset,
    y = marker_groups_subset
  )

  all_markers$group <- as.numeric(all_markers$group)

  return(list(
    markers = all_markers,
    mat = mat_subset,
    groups = marker_groups_subset
  ))
}

getMarkerNames <- function(data, all_markers) {
  all_markers$name <- data@misc$gene_annotations[all_markers$feature, "name"]
  all_markers <- all_markers |>
    dplyr::transmute(group = group, input = feature, name = name)

  rownames(all_markers) <- c()
  return(all_markers)
}

memoisedGetTopMarkerGenes <- memoise::memoise(
  getTopMarkerGenes,
  envir = .GlobalEnv,
  # cache_mem doesn't work because each request is run on a different process
  # so they don't share memory, so use cache_disk
  cache = cachem::cache_disk(
    dir="cache_marker_genes",
    destroy_on_finalize = FALSE
  ),
  # Ignore scdata changing (its size makes it a bad idea to hash) use cleanup_cache instead
  omit_args = c("data")
)

# Cleans up all the caches that depend on the seurat object
# should be run whenever the seurat object changes
cleanupMarkersCache <- function() {
  memoise::forget(memoisedGetTopMarkerGenes)
}
