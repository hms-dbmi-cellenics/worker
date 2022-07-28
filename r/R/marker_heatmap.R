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

  JSON_raw <- list()
  JSON_truncated <- list()

  JSON_raw$values <- mtx_res$rawExpression@x
  JSON_truncated$values <- mtx_res$truncatedExpression@x

  JSON_raw$index <- mtx_res$rawExpression@i
  JSON_truncated$index <- mtx_res$truncatedExpression@i

  JSON_raw$ptr <- mtx_res$rawExpression@p
  JSON_truncated$ptr <- mtx_res$truncatedExpression@p

  JSON_raw$size <- mtx_res$rawExpression@Dim
  JSON_truncated$size <- mtx_res$truncatedExpression@Dim


  res <- list()
  res$order <- names(stats)
  res$stats <- stats
  res$rawExpression <- JSON_raw
  res$truncatedExpression <- JSON_truncated

  return(res)
}
