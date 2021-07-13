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
  all_markers_original <- FindAllMarkers(data, min.pct = 0.25, min.diff.pct = 0.25)
  all_markers <- all_markers_original
  #Filtering out repeated genes to improve visualization, based on lowest p-value.
  #We could also use fold change.
  pvalueAggregate<-aggregate(all_markers$p_val_adj,by=list(all_markers$gene),FUN=min)
  minp <- pvalueAggregate$x[match(all_markers$gene,pvalueAggregate$Group.1)]
  all_markers$minp <- minp
  all_markers<-subset(all_markers,all_markers$p_val_adj<=all_markers$minp)
  all_markers <- all_markers[!duplicated(all_markers$gene), ]


  #Getting the indexes of the (nfeatures) features we want for each cluster.
  #This iteration will be as long as nclusters, so its acceptable.
  markerIndex <- rep(FALSE, length(all_markers$cluster))
  for (i in unique(all_markers$cluster)) {
    markerIndex <- markerIndex | all_markers$cluster == i
    markerIndex[which(all_markers$cluster == i, TRUE)[1 + nFeatures]:length(markerIndex)] <- FALSE
  }

  all_markers <- all_markers[markerIndex, ]

  df <- data@misc$gene_annotations
  genesSubset <- subset(df, toupper(df$input) %in% toupper(all_markers$gene))
  all_markers$name <- genesSubset[match(all_markers$gene,genesSubset$input),"name"]
  all_markers <- all_markers[,c("gene","name")]
  rownames(all_markers)<-c()
  colnames(all_markers)<-c("input","name")
  return(getExpressionValues(all_markers,data))
}
