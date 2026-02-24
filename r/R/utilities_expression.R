#' Extract expression values from Seurat object, add stats and format for UI
#'
#' @param data Seurat object
#' @param genes data.frame of genes of interest, with columns "input" and "name"
#' @param downsample_cell_ids vector. optional. If defined, only the expression is downsampled into these cells (every other cell is covered with 0's)
#'
#' @return list to send to the UI
#' @export
#'
getGeneExpression <- function(data, genes, downsample_cell_ids) {
  expression_values <- getExpressionValues(data, genes)

  ordered_gene_names <- ensure_is_list_in_json(colnames(expression_values))

  # getStats needs to use the real expression values (not downsampled) to extract correct stats
  stats <- getStats(expression_values)

  # If downsample_cell_ids exist, zero out cells not in the downsample set
  if (!missing(downsample_cell_ids)) {
    # Zero out all rows not in downsample_cell_ids
    all_cells <- data@meta.data$cells_id
    cells_to_zero <- setdiff(all_cells, downsample_cell_ids)
    zero_indices <- match(cells_to_zero, all_cells)
    expression_values[zero_indices, ] <- 0
  }

  # Format sparse matrix directly to JSON
  rawExpression <- toSparseJson(expression_values)

  return(list(
    orderedGeneNames = ordered_gene_names,
    stats = stats,
    rawExpression = rawExpression
  ))
}


#' Get raw gene expression values as sparse matrix
#'
#' @inheritParams getGeneExpression
#'
#' @return sparse matrix of expression values (cells x genes)
#' @export
#'
getExpressionValues <- function(data, genes) {
  mat <- data@assays$RNA$data

  # Subset to genes of interest and transpose (cells x genes)
  rawExpression <- Matrix::t(mat[unique(genes$input), , drop = FALSE])

  # Rename columns to display names
  symbol_idx <- match(colnames(rawExpression), genes$input)
  colnames(rawExpression) <- genes$name[symbol_idx]

  return(rawExpression)
}

#' Calculate the quantile truncation threshold for a vector
#'
#' @param x numeric vector
#' @param quantile_threshold numeric
#'
#' @return numeric truncation threshold
#'
getQuantileCap <- function(x, quantile_threshold) {
  lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
  i <- 0.01
  while (lim == 0 && i + quantile_threshold <= 1) {
    lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
    i <- i + 0.01
  }
  return(lim)
}

getStats <- function(expression) {
  stats_unsafe <- list(
    rawMean = unname(colMeans(expression, na.rm = TRUE)),
    rawStdev = unname(apply(expression, 2, sd, na.rm = TRUE)),
    truncatedMin = unname(apply(expression, 2, min, na.rm = TRUE)),
    truncatedMax = unname(apply(expression, 2, getQuantileCap, QUANTILE_THRESHOLD))
  )

  stats <- lapply(stats_unsafe, ensure_is_list_in_json)

  return(stats)
}

#' Extract sparse matrix attributes to mathJS-like sparse matrix format
#'
#' @param matrix sparse matrix
#'
#' @return list with sparse matrix attributes
#'
toSparseJson <- function(matrix) {
  response_unsafe <- list(
    values = matrix@x,
    index = matrix@i,
    ptr = matrix@p,
    size = matrix@Dim
  )

  response <- lapply(response_unsafe, ensure_is_list_in_json)

  return(response)
}
