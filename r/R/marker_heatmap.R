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
  data <- getClusters(req$body$type, req$body$config$resolution, data)

  markers <- getTopUpregulatedMarkers(data, req$body$nGenes)
  enids <- markers$feature

  getHeatmapExpression(data, enids)
}

# used by runMarkerHeatmap
getTopUpregulatedMarkers <- function(data, ntop) {

  markers <- presto::wilcoxauc(data, assay = "data", seurat_assay = "RNA")

  # filters out repeated genes picking lowest p-value
  markers <- markers %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(pval)) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(pval) %>%
    dplyr::slice_head(n = ntop) %>%
    dplyr::arrange(group)

  return(markers)
}




