#' Generates a marker heatmap
#'
#' @param req list with number of genes and cellsets in which to search for markers
#' @param data Seurat object
#'
#' @return list
#' @export
#'
runMarkerHeatmap <- function(req, data) {
  t_total_start <- Sys.time()
  message("runMarkerHeatmap: Starting")
  
  nFeatures <- req$body$nGenes
  cellSets <- req$body$cellSets$children
  cell_ids <- req$body$cellIds

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  message(sprintf("  → Calling memoisedGetTopMarkerGenes with %d features, %d cellsets, %d cell_ids", nFeatures, length(cellSets), length(cell_ids)))
  top_markers <- memoisedGetTopMarkerGenes(nFeatures, data, cell_sets_ids)
  
  t_marker_names_start <- Sys.time()
  message("  → Calling getMarkerNames")
  top_markers <- getMarkerNames(data, top_markers)
  message(sprintf("  ⏱️  After getMarkerNames: %.2fs", difftime(Sys.time(), t_marker_names_start, units = "secs")))

  t_gene_expr_start <- Sys.time()
  message("  → Calling getGeneExpression")
  geneExpression <- getGeneExpression(data, top_markers)
  message(sprintf("  ⏱️  After getGeneExpression: %.2fs", difftime(Sys.time(), t_gene_expr_start, units = "secs")))

  message(sprintf("✅ runMarkerHeatmap completed in %.2fs total", difftime(Sys.time(), t_total_start, units = "secs")))
  return(geneExpression)
}
