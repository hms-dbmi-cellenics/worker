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

  mtx_res <- list(
    rawExpression = sparsify(expression$rawExpression),
    truncatedExpression = sparsify(expression$truncatedExpression)
  )

  res <- list(
    order = names(stats),
    stats = stats,
    rawExpression = to_sparse_json(mtx_res$rawExpression),
    truncatedExpression = to_sparse_json(mtx_res$truncatedExpression)
  )

  return(res)
}


sparsify <- function(expression) {
  data.table::setnafill(expression, fill = 0)
  sparse_matrix <- Matrix::Matrix(Matrix::as.matrix(expression), sparse = T)

  return(sparse_matrix)

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
