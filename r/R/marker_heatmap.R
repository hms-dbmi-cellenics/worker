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

  top_markers <- getTopMarkerGenes(nFeatures, data, cellSets)
  top_markers <- getMarkerNames(data, top_markers)

  res <- getExpressionValues(top_markers, data)
  if (!length(res$rawExpression)) {
    stop(
      generateErrorMessage(
        ErrorCodes$NO_MARKER_GENES,
        "Couldn't get marker genes"
      )
    )
  }

  return(res)
}
