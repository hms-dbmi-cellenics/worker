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
  for(set in cellSets){
    filtered_cells <- intersect(set$cellIds,object_ids)
    data$custom[object_ids %in% filtered_cells] <- match(set$key,lapply(cellSets,function(x){x$key}))
  }

  all_markers <- presto::wilcoxauc(data, group_by="custom",assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)
  # Filtering out repeated genes to improve visualization, based on lowest p-value.
  # We could also use fold change.
  all_markers <- all_markers %>%
    group_by(feature) %>%
    slice(which.max(logFC))

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
