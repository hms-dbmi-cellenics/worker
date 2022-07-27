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
  res$rawExpression[is.na(res$rawExpression)] <- 0
  res$truncatedExpression[is.na(res$truncatedExpression)] <- 0
  mtx_res <- list()
  mtx_res$rawExpression <- Matrix::Matrix(Matrix::as.matrix(res$rawExpression),sparse=TRUE)
  mtx_res$truncatedExpression <- Matrix::Matrix(Matrix::as.matrix(res$truncatedExpression),sparse=TRUE)

  res <- formatExpression(res)
  return(res)
}
