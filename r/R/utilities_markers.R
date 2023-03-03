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
getTopMarkerGenes <- function(nFeatures, data, cellSets, cellSetsKeys = c(), aucMin = 0.3, pctInMin = 20, pctOutMax = 70) {
  data$marker_groups <- NA

  memoisedHola(1,2)

  object_ids <- data$cells_id
  for (i in seq_along(cellSets)) {
    set <- cellSets[[i]]
    filtered_cells <- intersect(set$cellIds, object_ids)
    data$marker_groups[object_ids %in% filtered_cells] <- i
  }

  all_markers <- presto::wilcoxauc(data,
    group_by = "marker_groups",
    assay = "data",
    seurat_assay = "RNA"
  )
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
