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

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  top_markers <- memoisedGetTopMarkerGenes(nFeatures, data, cell_sets_ids)
  top_markers <- getMarkerNames(data, top_markers)

  # Format response with only gene names
  ordered_gene_names <- ensure_is_list_in_json(unique(top_markers$name))
  response <- list(orderedGeneNames = ordered_gene_names)

  return(response)
}
