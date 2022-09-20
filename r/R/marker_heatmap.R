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
  mtx_res$rawExpression <- Matrix::Matrix(Matrix::as.matrix(expression$rawExpression),sparse=TRUE)
  mtx_res$truncatedExpression <- Matrix::Matrix(Matrix::as.matrix(expression$truncatedExpression),sparse=TRUE)

  res <- list(
    order = names(stats),
    stats = stats,
    rawExpression = to_sparse_json(mtx$rawExpression),
    truncatedExpression = to_sparse_json(mtx$truncatedExpression)
  )

  return(res)
}


to_sparse_json <- function(matrix) {
  sparse_json <-
    list(
      values = matrix@x,
      index = matrix@i,
      ptr = matrix@p,
      size = matrix@Dim
    )

  return(sparse_json)
}
