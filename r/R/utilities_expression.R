getGeneExpression <- function(data, genes) {

  expression_values <- getExpressionValues(genes, data)
  stats <- summaryStats(expression_values)

  mtx_res <- lapply(expression_values, sparsify)
  json_res <- lapply(mtx_res, toSparseJson)

  res <- list(
    order = as.list(names(stats)),
    stats = stats,
    rawExpression = json_res$rawExpression,
    truncatedExpression = json_res$truncatedExpression,
    zScore = json_res$zScore
  )

  return(res)
}


#' Returns expression values for selected genes
#'
#' @param genes - Must have names and input(ensmbl ids)
#' @param data
#'
#' @return
#' @export
#'
#' @examples
getExpressionValues <- function(genes, data) {
  quantile_threshold <- 0.95

  # Get the expression values for those genes in the corresponding matrix.
  rawExpression <- data.table::as.data.table(Matrix::t(data@assays$RNA@data[unique(genes$input), , drop = FALSE]))

  rawExpression[,cells_id := data@meta.data$cells_id]
  data.table::setorder(rawExpression, cols = "cells_id")

  # add back all filtered cell_ids as empty columns
  rawExpression <- rawExpression[
    data.table::CJ(cells_id = seq(0, max(data@meta.data$cells_id)), unique=TRUE),
    on=.(cells_id)
  ]

  rawExpression[, cells_id := NULL]

  symbol_idx <- match(colnames(rawExpression), genes$input)
  colnames(rawExpression) <- genes$name[symbol_idx]

  res <- list(
    rawExpression = rawExpression,
    truncatedExpression = truncateExpression(rawExpression, quantile_threshold),
    zScore = scaleExpression(rawExpression)
  )

  return(res)
}


truncateExpression <- function(rawExpression, quantile_threshold) {
  truncatedExpression <-
    rawExpression[, lapply(.SD, quantileTruncate, quantile_threshold), .SDcols = colnames(rawExpression)]

  return(truncatedExpression)
}


#' Truncate expression values after quantile threshold
#'
#' basically returns x > a => a else x. But takes into account the special case
#' when more than quantile_threshold percent of the data is 0. extending
#' the threshold bit by bit until it differs from 0. Therefore preserving some
#' dynamic range in the adjusted expression values.
#'
#' @param x numeric vector of raw expression values
#' @param quantile_threshold
#'
#' @return numeric vector of truncated expression values
#' @export
#'
quantileTruncate <- function(x, quantile_threshold) {
  lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
  i <- 0.01
  while (lim == 0 & i + quantile_threshold <= 1) {
    lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
    i <- i + 0.01
  }
  return(pmin(x, lim))
}


summaryStats <- function(data) {
  return(purrr::map2(data$rawExpression, data$truncatedExpression, summaryStatsAux))
}

summaryStatsAux <- function(raw, trunc) {
  return(list(
      rawMean = mean(raw, na.rm = TRUE),
      rawStdev = sd(raw, na.rm = TRUE),
      truncatedMin = min(trunc, na.rm = TRUE),
      truncatedMax = max(trunc, na.rm = TRUE)
  ))
}


scaleExpression <- function(rawExpression) {
  scaledExpression <- rawExpression[, lapply(.SD, scale), .SDcols = colnames(rawExpression)]
  return(scaledExpression)
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
