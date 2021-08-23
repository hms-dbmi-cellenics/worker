#' Generates a marker heatmap
#'
#' @param req
#' @param data
#'
#' @return
#' @export
#'
#' @examples
runMarkerHeatmap <- function(req, data) {
  nFeatures <- req$body$nGenes
  cellSets <- req$body$cellSets$children
  data$custom <- NA

  object_ids <- data$cells_id
  for (i in 1:length(cellSets)){
    set <- cellSets[[i]]
    filtered_cells <- intersect(set$cellIds, object_ids)
    data$custom[object_ids %in% filtered_cells] <- i
  }

  all_markers <- presto::wilcoxauc(data, group_by="custom", assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)
  # Filtering out repeated genes to avoid displaying the same genes for two groups, based on lowest p-value
  all_markers <- all_markers %>%
    dplyr::filter(logFC > 0) %>%
    group_by(feature) %>%
    slice(which.min(pval))

  all_markers <- all_markers %>%
    group_by(group) %>%
    arrange(pval) %>%
    dplyr::slice_head(n = nFeatures) %>%
    arrange(group)

  df <- data@misc$gene_annotations
  genesSubset <- subset(df, toupper(df$input) %in% toupper(all_markers$feature))
  all_markers$name <- genesSubset[match(all_markers$feature, genesSubset$input), "name"]
  all_markers <- all_markers[, c("feature", "name")]
  rownames(all_markers) <- c()
  colnames(all_markers) <- c("input", "name")
  return(getExpressionValues(all_markers, data))
}
