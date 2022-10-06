#' Generates a marker heatmap
#'
#' @param req list with number of genes and cellsets in which to search for markers
#' @param data Seurat object
#'
#' @return list
#' @export
#'
runMarkerHeatmap <- function(req, data) {
  nFeatures <- req$body$nGenes
  cellSets <- req$body$cellSets$children

  top_markers <- getTopMarkerGenes(nFeatures, data, cellSets)
  top_markers <- getMarkerNames(data, top_markers)

  return(getGeneExpression(data, top_markers))
}
