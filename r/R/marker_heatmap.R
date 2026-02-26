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

  # Format response to match getGeneExpression structure, but with only markers
  message("  → Formatting response")
  ordered_gene_names <- ensure_is_list_in_json(unique(top_markers$name))
  response <- list(
    orderedGeneNames = ordered_gene_names
  )
  message(sprintf("✅ runMarkerHeatmap completed in %.2fs total", difftime(Sys.time(), t_total_start, units = "secs")))
  return(response)
}
