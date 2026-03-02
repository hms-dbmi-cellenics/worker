#' Expand sparse matrix to full cell dimensions using cell IDs
#'
#' Converts sparse matrix to triplet format and remaps row indices to actual
#' cell IDs, creating a new sparse matrix with dimensions covering all cells.
#'
#' @param expression_values sparse matrix (subset of cells x genes)
#' @param all_cell_ids vector of cell IDs (0-based) corresponding to rows in expression_values
#'
#' @return sparse matrix with full dimensions (n_full_cells x genes)
#'
expandMatrixToCellIDs <- function(expression_values, all_cell_ids) {
  # Convert to triplet format
  trip <- as(expression_values, "TsparseMatrix")
  
  # Map row indices to their actual cell_ids
  # trip@i is 0-based row position in current matrix, all_cell_ids gives the cell_id for that row
  new_i <- all_cell_ids[trip@i + 1] + 1  # cell_ids are 0-based, convert to 1-based for sparseMatrix
  new_j <- trip@j + 1  # Convert to 1-based for sparseMatrix
  
  # Calculate full dimensions based on cell_ids
  max_cell_id <- max(all_cell_ids)
  n_full_cells <- max_cell_id + 1
  n_genes <- ncol(expression_values)
  
  # Create new sparse matrix with full dimensions (all cells)
  expanded_matrix <- Matrix::sparseMatrix(
    i = new_i,
    j = new_j,
    x = trip@x,
    dims = c(n_full_cells, n_genes),
    giveCsparse = TRUE
  )
  
  return(expanded_matrix)
}

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

  ordered_gene_names <- ensure_is_list_in_json(colnames(expression_values))

  # getStats uses the expression values for all cells
  stats <- getStats(expression_values)

  # Expand matrix to full dimensions based on cell_ids
  # This ensures row indices correspond to cell_id positions in the full dataset
  all_cell_ids <- data@meta.data$cells_id
  expression_values <- expandMatrixToCellIDs(expression_values, all_cell_ids)

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

#' Calculate the quantile truncation threshold for sparse matrix columns
#'
#' @param x sparse matrix
#' @param quantile_threshold numeric
#'
#' @return numeric vector of truncation thresholds, one per column
#'
getQuantileCap <- function(x, quantile_threshold) {
  
  # Calculate quantiles for all columns at once
  lims <- sparseMatrixStats::colQuantiles(x, probs = quantile_threshold, na.rm = TRUE, drop = TRUE)
  
  # Find columns where quantile is 0 and needs adjustment
  zero_cols <- which(lims == 0)
  
  if (length(zero_cols) > 0) {
    # Iterate through quantile thresholds for zero columns only
    q_threshold <- quantile_threshold + 0.01
    
    while (q_threshold <= 1 && length(zero_cols) > 0) {
      # Calculate new quantiles only for columns that are still 0
      new_lims <- sparseMatrixStats::colQuantiles(
        x[, zero_cols, drop = FALSE], 
        probs = q_threshold, 
        na.rm = TRUE, 
        drop = TRUE)
      
      # Update the limits
      lims[zero_cols] <- new_lims
      
      # Update which columns are still 0
      zero_cols <- zero_cols[new_lims == 0]
      
      q_threshold <- q_threshold + 0.01
    }
  }
  
  return(as.numeric(lims))
}

getStats <- function(expression) {
  # Matrix::colMeans works directly with sparse matrices efficiently
  mean_vals <- unname(Matrix::colMeans(expression, na.rm = TRUE))
  
  # Optimize stdev calculation for sparse matrices
  # Formula: sqrt(sum((x - mean)^2) / (n - 1))
  # Vectorized approach is much faster than apply(sd)
  n <- nrow(expression)
  # Center the matrix by subtracting column means
  centered <- expression - Matrix::Matrix(rep(mean_vals, n), nrow = n, byrow = TRUE)
  # Square the centered values
  centered_sq <- centered^2
  # Sum of squares per column
  sum_sq <- Matrix::colSums(centered_sq, na.rm = TRUE)
  # Standard deviation = sqrt(sum_sq / (n - 1))
  stdev_vals <- unname(sqrt(sum_sq / (n - 1)))
  
  min_vals <- unname(apply(expression, 2, min, na.rm = TRUE))
  
  max_vals <- getQuantileCap(expression, QUANTILE_THRESHOLD)
  
  stats_unsafe <- list(
    rawMean = mean_vals,
    rawStdev = stdev_vals,
    truncatedMin = min_vals,
    truncatedMax = max_vals
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
