#' getTopMarkerGenes
#'
#' @param nFeatures
#' @param data
#' @param cellSets
#' @param aucMin
#' @param pctInMin
#' @param pctOutMax
#'
#' @return Returns top marker genes as calculated by presto
#' @export
#'
#' @examples
getTopMarkerGenes <- function(nFeatures, data, cellSets, aucMin = 0.3, pctInMin = 20, pctOutMax = 70) {
  data$marker_groups <- NA

  object_ids <- data$cells_id
  for (i in seq_along(cellSets)) {
    set <- cellSets[[i]]
    filtered_cells <- intersect(set$cellIds, object_ids)
    data$marker_groups[object_ids %in% filtered_cells] <- i
  }

  all_markers <- presto::wilcoxauc(data, group_by = "marker_groups", assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)

  # may not return nFeatures markers per cluster if values are too stringent
  filtered_markers <- all_markers %>%
    dplyr::filter(logFC > 0 &
      auc >= aucMin &
      pct_in >= pctInMin &
      pct_out <= pctOutMax) %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(pval))

  top_markers <- filtered_markers %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::arrange(group)

  message(sprintf("%d markers selected", nrow(top_markers)))
  return(top_markers)
}

getMarkerNames <- function(data, all_markers) {
  all_markers$name <- data@misc$gene_annotations[all_markers$feature, "name"]
  all_markers <- all_markers %>% dplyr::transmute(input = feature, name = name)
  rownames(all_markers) <- c()
  return(all_markers)
}
