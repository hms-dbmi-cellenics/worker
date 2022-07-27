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

  expression <- getExpressionValues(top_markers, data)
  stats <- formatExpressionMtx(expression)

  expression$rawExpression[is.na(expression$rawExpression)] <- 0
  expression$truncatedExpression[is.na(expression$truncatedExpression)] <- 0

  mtx_res <- list()
  mtx_res$rawExpression <- Matrix::Matrix(Matrix::as.matrix(t(expression$rawExpression)),sparse=TRUE)
  mtx_res$truncatedExpression <- Matrix::Matrix(Matrix::as.matrix(t(expression$truncatedExpression)),sparse=TRUE)

  res <- list()
  res$order <- names(stats)
  res$stats <- stats
  res$rawExpression <- mtx_res$rawExpression
  res$truncatedExpression <- mtx_res$truncatedExpression

  return(res)
}
