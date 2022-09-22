#' Extract expression values from Seurat object, add stats and format for UI
#'
#' @param data Seurat object
#' @param genes character vector og genes to extract
#'
#' @return list
#' @export
#'
getGeneExpression <- function(data, genes) {
  expression_values <- getExpressionValues(data, genes)
  stats <- summaryStats(expression_values)

  mtx_res <- lapply(expression_values, sparsify)
  json_res <- lapply(mtx_res, toSparseJson)

  res <- list(
    order = names(stats),
    stats = stats,
    rawExpression = json_res$rawExpression,
    truncatedExpression = json_res$truncatedExpression,
    zScore = json_res$zScore
  )

  return(res)
}


#' Get raw, truncated and scaled gene expression values
#'
#' @inheritParams getGeneExpression
#'
#' @return list of expression values
#' @export
#'
getExpressionValues <- function(data, genes) {

  rawExpression <- getRawExpression(data, genes)

  res <- list(
    rawExpression = rawExpression,
    truncatedExpression = truncateExpression(rawExpression),
    zScore = scaleExpression(rawExpression)
  )

  return(res)
}

#' Extract raw expression values for a list of genes
#'
#' The expression matrix is transposed to accommodate CSC sparse matrix format
#' used in the UI. Filtered cells are added as empty rows.
#'
#' @inheritParams getGeneExpression
#'
#' @return data.table of raw expression values
#' @export
#'
getRawExpression <- function(data, genes) {
  rawExpression <-
    Matrix::t(data@assays$RNA@data[unique(genes$input), , drop = FALSE])
  rawExpression <-   data.table::as.data.table(rawExpression)

  rawExpression <- completeRawExpression(rawExpression, data@meta.data$cells_id)

  symbol_idx <- match(colnames(rawExpression), genes$input)
  colnames(rawExpression) <- genes$name[symbol_idx]

  return(rawExpression)

}

#' Adds an empty row for every filtered cell
#'
#' The UI infers cell_id by the index of the cell in the matrix, which means that
#' filtered cells have to be added back to the table, as empty rows. When converted
#' to sparse format, they do not take up space.
#'
#' @param rawExpression data.table
#' @param cell_ids integer vector
#'
#' @return complete raw expression data.table
#' @export
#'
completeRawExpression <- function(rawExpression, cell_ids) {
  rawExpression[, cell_ids := cell_ids]
  data.table::setorder(rawExpression, cols = "cell_ids")

  # add back all filtered cell_ids as empty columns
  rawExpression <-
    rawExpression[data.table::CJ(cell_ids = seq(0, max(cell_ids)),
                                 unique = TRUE),
                  on = .(cell_ids)]

  rawExpression[, cell_ids := NULL]

  return(rawExpression)
}


#' Truncates expression values for all genes in data.table
#'
#' @param rawExpression data.table
#'
#' @return data.table of truncated gene expression values
#' @export
#'
truncateExpression <- function(rawExpression) {
  truncatedExpression <-
    rawExpression[, lapply(.SD, quantileTruncate, QUANTILE_THRESHOLD), .SDcols = colnames(rawExpression)]

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
#' @param quantile_threshold numeric
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

#' Calculate z-score of gene expression values
#'
#' Centers values to the mean and scales to the standard deviation.
#'
#' @param rawExpression data.table
#'
#' @return
#' @export
#'
scaleExpression <- function(rawExpression) {
  scaledExpression <-
    rawExpression[, lapply(.SD, scale), .SDcols = colnames(rawExpression)]
  return(scaledExpression)
}


summaryStats <- function(data) {
  return(purrr::map2(
    data$rawExpression,
    data$truncatedExpression,
    summaryStatsAux
  ))
}

summaryStatsAux <- function(raw, trunc) {
  return(list(
    rawMean = mean(raw, na.rm = TRUE),
    rawStdev = sd(raw, na.rm = TRUE),
    truncatedMin = min(trunc, na.rm = TRUE),
    truncatedMax = max(trunc, na.rm = TRUE)
  ))
}


#' Convert data.table to CSC sparse matrix
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
  sparse_matrix <-
    Matrix::Matrix(Matrix::as.matrix(expression), sparse = T)

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
