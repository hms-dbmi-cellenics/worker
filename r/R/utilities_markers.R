#' getTopMarkerGenes
#'
#' Uses presto::wilcoxauc to find the marker genes that distinguish the
#' cellsets. It then filters the list of genes up to nFeatures using reasonable
#' defaults.
#'
#' @param nFeatures int number of marker genes to get
#' @param data SeuratObject
#' @param cellSets list of cellsets to split for marker gene selection
#' @param aucMin min area under the wilcoxon test's ROC for a gene to be
#'  considered a marker
#' @param pctInMin min percentage of cells in cellset that have to express a
#'  gene for it to be considered a marker
#' @param pctOutMax max percentage of cells outside cellset that can express a
#'  gene for it to be considered a marker
#'
#' @return data.frame of top marker genes
#' @export
#'
getTopMarkerGenes <- function(nFeatures, data, cellSetsIds, aucMin = 0.3, pctInMin = 20, pctOutMax = 70) {
  t_start <- Sys.time()
  message("Running getTopMarkerGenes")
  
  object_ids <- data$cells_id
  marker_groups <- rep(NA, length(object_ids))

  for (i in seq_along(cellSetsIds)) {
    filtered_cells <- intersect(cellSetsIds[[i]], object_ids)
    marker_groups[object_ids %in% filtered_cells] <- i
  }
  message(sprintf("  ⏱️  Assigned marker groups: %.2fs", difftime(Sys.time(), t_start, units = "secs")))

  # for speed: take at most 1000 cells per cluster
  set.seed(0)
  t_sample_start <- Sys.time()
  
  # Create temp data frame for sampling
  temp_meta <- data@meta.data |>
    dplyr::mutate(marker_groups = marker_groups)
  
  keep_cell_ids <- temp_meta |>
    dplyr::group_by(marker_groups) |>
    dplyr::slice_sample(n = 1000, replace = FALSE) |>
    dplyr::ungroup() |>
    dplyr::pull(cells_id)
  message(sprintf("  ⏱️  Sampling cells (%d cells selected): %.2fs", length(keep_cell_ids), difftime(Sys.time(), t_sample_start, units = "secs")))
  
  # Extract and subset matrix
  t_matrix_start <- Sys.time()
  mat <- data@assays$RNA$data
  keep_indices <- match(keep_cell_ids, object_ids)
  mat_subset <- mat[, keep_indices]
  marker_groups_subset <- marker_groups[keep_indices]
  message(sprintf("  ⏱️  Matrix subsetting: %.2fs (subset dim: %d x %d)", difftime(Sys.time(), t_matrix_start, units = "secs"), nrow(mat_subset), ncol(mat_subset)))

  t_wilcox_start <- Sys.time()
  all_markers <- presto::wilcoxauc(
    mat_subset,
    y = marker_groups_subset
  )
  message(sprintf("  ⏱️  Wilcoxauc test: %.2fs (%d results)", difftime(Sys.time(), t_wilcox_start, units = "secs"), nrow(all_markers)))
  
  all_markers$group <- as.numeric(all_markers$group)

  # may not return nFeatures markers per cluster if values are too stringent
  t_filter_start <- Sys.time()
  filtered_markers <- all_markers %>%
    dplyr::filter(logFC > 0 &
      auc >= aucMin &
      pct_in >= pctInMin &
      pct_out <= pctOutMax) %>%
    dplyr::arrange(pval) %>%
    dplyr::distinct(feature, .keep_all = TRUE) %>%
    dplyr::arrange(-logFC)
  message(sprintf("  ⏱️  Filtering markers: %.2fs (%d markers after filtering)", difftime(Sys.time(), t_filter_start, units = "secs"), nrow(filtered_markers)))

  # Get top nFeatures per group without additional sorting
  t_top_start <- Sys.time()
  top_markers <- filtered_markers %>%
    dplyr::group_by(group) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(group)
  message(sprintf("  ⏱️  Selecting top markers: %.2fs (%d final markers)", difftime(Sys.time(), t_top_start, units = "secs"), nrow(top_markers)))
  
  message(sprintf("✅ %d markers selected in %.2fs total", nrow(top_markers), difftime(Sys.time(), t_start, units = "secs")))
  return(top_markers)
}

getMarkerNames <- function(data, all_markers) {
  t_start <- Sys.time()
  
  all_markers$name <- data@misc$gene_annotations[all_markers$feature, "name"]
  t_transmute_start <- Sys.time()
  all_markers <- all_markers %>% dplyr::transmute(group = group, input = feature, name = name)
  message(sprintf("  ⏱️  getMarkerNames: %.2fs", difftime(Sys.time(), t_transmute_start, units = "secs")))
  
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
