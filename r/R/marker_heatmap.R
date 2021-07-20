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
  data <- getClusters(req$body$type, req$body$config$resolution, data)

  all_markers <- presto::wilcoxauc(data, assay = "data", seurat_assay = "RNA")
  # Filtering out repeated genes to improve visualization, based on lowest p-value.
  # We could also use fold change.
 all_markers <- all_markers %>%  
 group_by(feature) %>% 
 slice(which.min(pval))

  nFeatures <- as.integer(30 / (as.integer(max(all_markers$group)) + 1))
  all_markers <- all_markers %>%
   group_by(group) %>%
    arrange(padj) %>%
    slice_head(n=5)

  df <- data@misc$gene_annotations
  genesSubset <- subset(df, toupper(df$input) %in% toupper(all_markers$feature))
  all_markers$name <- genesSubset[match(all_markers$feature, genesSubset$input), "name"]
  all_markers <- all_markers[, c("feature", "name")]
  rownames(all_markers) <- c()
  colnames(all_markers) <- c("input", "name")
  return(getExpressionValues(all_markers, data))
}
