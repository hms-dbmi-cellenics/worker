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
  message("\n\n*** CORES")
  message(data.table::getDTthreads())
  message("*** CORES\n\n")

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

  adjGeneExpression <- truncateExpression(rawExpression, quantile_threshold)

  return(list(rawExpression = rawExpression, truncatedExpression = adjGeneExpression))
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


#' formatExpression
#'
#' @param data
#'
#' @return Formats expression values for the UI.
#' @export
#'
#' @examples
formatExpression <- function(data) {
  return(purrr::map2(data$rawExpression, data$truncatedExpression, formatExpressionAux))
}

formatExpressionAux <- function(raw, trunc) {
  return(list(
    rawExpression = list(
      mean = mean(raw, na.rm = TRUE),
      stdev = sd(raw, na.rm = TRUE),
      expression = raw
    ),
    truncatedExpression = list(
      min = min(trunc, na.rm = TRUE),
      max = max(trunc, na.rm = TRUE),
      expression = trunc
    )
  ))
}

formatExpressionMtx <- function(data) {
  return(purrr::map2(data$rawExpression, data$truncatedExpression, formatExpressionMtxAux))
}

formatExpressionMtxAux <- function(raw, trunc) {
  return(list(
      rawMean = mean(raw, na.rm = TRUE),
      rawStdev = sd(raw, na.rm = TRUE),
      truncatedMin = min(trunc, na.rm = TRUE),
      truncatedMax = max(trunc, na.rm = TRUE)
  ))
}


