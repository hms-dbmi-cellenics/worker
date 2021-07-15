#' Generates a marker heatmap
#'
#' @param req
#' @param data
#'
#' @return
#' @export
#'
#' @examples
runMarkerHeatmap <- function(req,data) {
  nFeatures <- req$body$nGene
  data <- getClusters(req$body$type,req$body$config$resolution,data)
  all_markers<-presto::wilcoxauc(data,assay = "data", seurat_assay = "RNA")

  #Filtering out repeated genes to improve visualization, based on lowest p-value.
  #We could also use fold change.
  pvalueAggregate<-aggregate(all_markers$padj,by=list(all_markers$feature),FUN=min)
  minp <- pvalueAggregate$x[match(all_markers$feature,pvalueAggregate$Group.1)]
  all_markers$minp <- minp
  all_markers<-subset(all_markers,all_markers$padj<=all_markers$minp)
  all_markers <- all_markers[!duplicated(all_markers$feature), ]

  all_markers <- all_markers %>%  group_by(group) %>% arrange(padj) %>% arrange(group) %>% dplyr::filter(row_number() %in% c(1:nFeatures))

  df <- data@misc$gene_annotations
  genesSubset <- subset(df, toupper(df$input) %in% toupper(all_markers$feature))
  all_markers$name <- genesSubset[match(all_markers$feature,genesSubset$input),"name"]
  all_markers <- all_markers[,c("feature","name")]
  rownames(all_markers)<-c()
  colnames(all_markers)<-c("input","name")
  return(getExpressionValues(all_markers,data))
}



