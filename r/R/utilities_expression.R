#' Extract expression values from Seurat object, add stats and format for UI
#'
#' @param data Seurat object
#' @param genes data.frame of genes of interest, with columns "input" and "name"
#'
#' @return list to send to the UI
#' @export
#'
getGeneExpression <- function(data, genes) {
  expression_values <- getExpressionValues(data, genes)
  stats <- summaryStats(expression_values)

  expression_values <-
    lapply(expression_values, formatExpression, data@meta.data$cells_id)

  return(list(
    order = as.list(names(stats)),
    stats = stats,
    rawExpression = expression_values$rawExpression,
    truncatedExpression = expression_values$truncatedExpression,
    zScore = expression_values$zScore
  ))
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

  return(list(
    rawExpression = rawExpression,
    truncatedExpression = truncateExpression(rawExpression),
    zScore = scaleExpression(rawExpression)
  ))
}

#' Extract raw expression values for a list of genes
#'
#' The expression matrix is transposed to accommodate the CSC sparse matrix
#' format used in the UI.
#'
#' @inheritParams getGeneExpression
#'
#' @return data.table of raw expression values
#' @export
#'
getRawExpression <- function(data, genes) {
  rawExpression <-
    Matrix::t(data@assays$RNA@data[unique(genes$input), , drop = FALSE])
  rawExpression <- data.table::as.data.table(rawExpression)

  symbol_idx <- match(colnames(rawExpression), genes$input)
  colnames(rawExpression) <- genes$name[symbol_idx]

  return(rawExpression)
}

#' Adds an empty row for every filtered cell
#'
#' The UI infers cell_id by the index of the cell in the matrix, which means
#' that filtered cells have to be added back to the table as empty rows. When
#' converted to sparse format, they do not take up space. We use the max of the
#' cell_ids that were filtered, which means this table will not contain cells
#' above that. But it does not change the index of the cells below, and the UI
#' does not care.
#'
#' @param expression data.table
#' @param cell_ids integer vector
#'
#' @return row complete expression data.table
#' @export
#'
completeExpression <- function(expression, cell_ids) {
  expression[, cell_ids := cell_ids]
  data.table::setorder(expression, cols = "cell_ids")

  # add back all filtered cells as empty rows.
  expression <-
    expression[data.table::CJ(
      cell_ids = seq(0, max(cell_ids)),
      unique = TRUE
    ),
    on = .(cell_ids)
    ]

  expression[, cell_ids := NULL]

  return(expression)
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
    rawExpression[, lapply(.SD, quantileTruncate, QUANTILE_THRESHOLD),
      .SDcols = colnames(rawExpression)
    ]

  return(truncatedExpression)
}


#' Truncate expression values after quantile threshold
#'
#' Basically returns x > a => a else x. But takes into account the special case
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
  while (lim == 0 && i + quantile_threshold <= 1) {
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
#' @return data.table of scaled expression values
#' @export
#'
scaleExpression <- function(rawExpression) {
  # scale returns a matrix. need to convert to vector avoid changing colnames
  scaledExpression <-
    rawExpression[, lapply(.SD, \(x) as.vector(scale(x))),
      .SDcols = colnames(rawExpression)
    ]
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
  return(list(
    values = matrix@x,
    index = matrix@i,
    ptr = matrix@p,
    size = matrix@Dim
  ))
}


#' Format expression data.table as mathJS json
#'
#' @param expression data.table with expression values
#' @param cell_ids int vector with filtered cell_ids
#'
#' @return list of formatted expression table
#' @export
#'
formatExpression <- function(expression, cell_ids) {
  expression <- completeExpression(expression, cell_ids)
  expression <- sparsify(expression)
  expression <- toSparseJson(expression)

  return(expression)
}
