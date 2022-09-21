#' Generates a marker heatmap
#'
#' @param req
#' @param data
#'
#' @return list
#' @export
#'
runMarkerHeatmap <- function(req, data) {
  nFeatures <- req$body$nGenes
  cellSets <- req$body$cellSets$children

  top_markers <- getTopMarkerGenes(nFeatures, data, cellSets)
  top_markers <- getMarkerNames(data, top_markers)

  return(getGeneExpression(data, top_markers))
}


getGeneExpression <- function(data, genes) {

  expression_values <- getExpressionValues(genes, data)
  stats <- summaryStats(expression_values)

  mtx_res <- list(
    rawExpression = sparsify(expression_values$rawExpression),
    truncatedExpression = sparsify(expression_values$truncatedExpression)
  )

  res <- list(
    order = names(stats),
    stats = stats,
    rawExpression = toSparseJson(mtx_res$rawExpression),
    truncatedExpression = toSparseJson(mtx_res$truncatedExpression)
  )

  return(res)
}


#' convert data.table to CSC sparse matrix
#'
#' NAs are replaced by zeroes by reference. Then coerced to sparse matrix.
#'
#' @param expression data.table
#'
#' @return dgCmatrix
#' @export
#'
sparsify <- function(expression) {
  data.table::setnafill(expression, fill = 0)
  sparse_matrix <- Matrix::Matrix(Matrix::as.matrix(expression), sparse = T)

  return(sparse_matrix)
}

#' extract sparse matrix attributes to mathJS-like sparse matrix format
#'
#' @param matrix
#'
#' @return list with sparse matrix attributes
#' @export
#'
toSparseJson <- function(matrix) {
  sparse_json <-
    list(
      values = matrix@x,
      index = matrix@i,
      ptr = matrix@p,
      size = matrix@Dim
    )

  return(sparse_json)
}
