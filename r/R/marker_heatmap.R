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
  cell_ids <- req$body$cellIds

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  top_markers <- memoisedGetTopMarkerGenes(nFeatures, data, cellSets, cell_sets_ids)
  top_markers <- getMarkerNames(data, top_markers)

  geneExpression <- getGeneExpression(data, top_markers, cell_ids)

  return(geneExpression)
}
