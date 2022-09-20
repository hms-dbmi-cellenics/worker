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
  library(data.table)
  quantile_threshold <- 0.95

  # Get the expression values for those genes in the corresponding matrix.
  geneExpression <- data.table::as.data.table(Matrix::t(data@assays$RNA@data[unique(genes$input), , drop = FALSE]))

  geneExpression[,cells_id := data@meta.data$cells_id]
  data.table::setorder(geneExpression, cols = "cells_id")

  # add back all filtered cell_ids
  geneExpression <- geneExpression[
    data.table::CJ(cells_id = seq(0, max(data@meta.data$cells_id)), unique=TRUE),
    on=.(cells_id)
  ]

  geneExpression[, cells_id := NULL]

  symbol_idx <- match(colnames(geneExpression), genes$input)
  colnames(geneExpression) <- genes$name[symbol_idx]

  adjGeneExpression <- truncateExpression(geneExpression)

  return(list(rawExpression = geneExpression, truncatedExpression = adjGeneExpression))
}


truncateExpression <- function(geneExpression) {
  adj_expresion <-
    geneExpression[, lapply(.SD, q_transform), .SDcols = colnames(geneExpression)]

  return(adj_expresion)
}


q_transform <- function(x) {
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


