#' Extract expression values from Seurat object, add stats and format for UI
#'
#' @param data Seurat object
#' @param genes data.frame of genes of interest, with columns "input" and "name"
#'
#' @return list to send to the UI
#' @export
#'
getGeneExpression <- function(data, genes) {
  t_start <- Sys.time()
  message("getGeneExpression: Starting")
  
  t_expr_start <- Sys.time()
  expression_values <- getExpressionValues(data, genes)
  message(sprintf("  ⏱️  getExpressionValues: %.2fs (dim: %d x %d)", difftime(Sys.time(), t_expr_start, units = "secs"), nrow(expression_values), ncol(expression_values)))

  t_names_start <- Sys.time()
  ordered_gene_names <- ensure_is_list_in_json(colnames(expression_values))
  message(sprintf("  ⏱️  Ordered gene names: %.2fs", difftime(Sys.time(), t_names_start, units = "secs")))

  # getStats uses the expression values for all cells
  t_stats_start <- Sys.time()
  stats <- getStats(expression_values)
  message(sprintf("  ⏱️  getStats: %.2fs", difftime(Sys.time(), t_stats_start, units = "secs")))

  # Expand matrix to full dimensions based on cell_ids
  # This ensures row indices correspond to cell_id positions in the full dataset
  t_expand_start <- Sys.time()
  all_cell_ids <- data@meta.data$cells_id
  max_cell_id <- max(all_cell_ids)
  n_full_cells <- max_cell_id + 1
  n_genes <- ncol(expression_values)
  
  # Convert to triplet format
  trip <- as(expression_values, "TsparseMatrix")
  
  # Map row indices to their actual cell_ids
  # trip@i is 0-based row position in current matrix, all_cell_ids gives the cell_id for that row
  new_i <- all_cell_ids[trip@i + 1] + 1  # cell_ids are 0-based, convert to 1-based for sparseMatrix
  new_j <- trip@j + 1  # Convert to 1-based for sparseMatrix
  
  # Create new sparse matrix with full dimensions (all cells, no downsampling)
  expression_values <- Matrix::sparseMatrix(
    i = new_i,
    j = new_j,
    x = trip@x,
    dims = c(n_full_cells, n_genes),
    giveCsparse = TRUE
  )
  message(sprintf("  ⏱️  Matrix expansion: %.2fs (now %d x %d)", difftime(Sys.time(), t_expand_start, units = "secs"), nrow(expression_values), ncol(expression_values)))

  # Format sparse matrix directly to JSON
  t_json_start <- Sys.time()
  message("  ⚠️  Starting toSparseJson conversion...")
  rawExpression <- toSparseJson(expression_values)
  message(sprintf("  ⏱️  toSparseJson: %.2fs", difftime(Sys.time(), t_json_start, units = "secs")))

  message(sprintf("✅ getGeneExpression completed in %.2fs total", difftime(Sys.time(), t_start, units = "secs")))
  
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
  t_start <- Sys.time()
  
  t_matrix_start <- Sys.time()
  mat <- data@assays$RNA$data
  message(sprintf("  📊 Matrix loaded: %.2fs (dim: %d x %d)", difftime(Sys.time(), t_matrix_start, units = "secs"), nrow(mat), ncol(mat)))

  # Subset to genes of interest and transpose (cells x genes)
  t_subset_start <- Sys.time()
  rawExpression <- Matrix::t(mat[unique(genes$input), , drop = FALSE])
  message(sprintf("  ✂️  Matrix subset & transpose: %.2fs (result: %d x %d)", difftime(Sys.time(), t_subset_start, units = "secs"), nrow(rawExpression), ncol(rawExpression)))

  # Rename columns to display names
  t_rename_start <- Sys.time()
  symbol_idx <- match(colnames(rawExpression), genes$input)
  colnames(rawExpression) <- genes$name[symbol_idx]
  message(sprintf("  🏷️  Column rename: %.2fs", difftime(Sys.time(), t_rename_start, units = "secs")))

  message(sprintf("  → getExpressionValues total: %.2fs", difftime(Sys.time(), t_start, units = "secs")))
  return(rawExpression)
}

# Fixed vectorized version
getQuantileCap_vectorized <- function(x, quantile_threshold) {
  
  # Calculate quantiles for all columns at once
  lims <- sparseMatrixStats::colQuantiles(x, probs = quantile_threshold, na.rm = TRUE, drop = TRUE)
  
  # Find columns where quantile is 0 and needs adjustment
  zero_cols <- which(lims == 0)
  
  if (length(zero_cols) > 0) {
    # Iterate through quantile thresholds for zero columns only
    q_threshold <- quantile_threshold + 0.01
    
    while (q_threshold <= 1 && length(zero_cols) > 0) {
      # Calculate new quantiles only for columns that are still 0
      new_lims <- sparseMatrixStats::colQuantiles(x[, zero_cols, drop = FALSE], 
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
  t_start <- Sys.time()
  
  # Matrix::colMeans works directly with sparse matrices efficiently
  t_mean_start <- Sys.time()
  mean_vals <- unname(Matrix::colMeans(expression, na.rm = TRUE))
  message(sprintf("  📈 colMeans: %.2fs", difftime(Sys.time(), t_mean_start, units = "secs")))
  
  # Optimize stdev calculation for sparse matrices
  # Formula: sqrt(sum((x - mean)^2) / (n - 1))
  # Vectorized approach is much faster than apply(sd)
  t_stdev_start <- Sys.time()
  n <- nrow(expression)
  # Center the matrix by subtracting column means
  centered <- expression - Matrix::Matrix(rep(mean_vals, n), nrow = n, byrow = TRUE)
  # Square the centered values
  centered_sq <- centered^2
  # Sum of squares per column
  sum_sq <- Matrix::colSums(centered_sq, na.rm = TRUE)
  # Standard deviation = sqrt(sum_sq / (n - 1))
  stdev_vals <- unname(sqrt(sum_sq / (n - 1)))
  message(sprintf("  📊 stdev (vectorized): %.2fs", difftime(Sys.time(), t_stdev_start, units = "secs")))
  
  t_min_start <- Sys.time()
  min_vals <- unname(apply(expression, 2, min, na.rm = TRUE))
  message(sprintf("  📉 min: %.2fs", difftime(Sys.time(), t_min_start, units = "secs")))
  
  t_max_start <- Sys.time()
  max_vals <- getQuantileCap_vectorized(expression, QUANTILE_THRESHOLD)
  message(sprintf("  📈 max (quantile): %.2fs", difftime(Sys.time(), t_max_start, units = "secs")))
  
  stats_unsafe <- list(
    rawMean = mean_vals,
    rawStdev = stdev_vals,
    truncatedMin = min_vals,
    truncatedMax = max_vals
  )

  t_json_start <- Sys.time()
  stats <- lapply(stats_unsafe, ensure_is_list_in_json)
  message(sprintf("  🔄 JSON formatting: %.2fs", difftime(Sys.time(), t_json_start, units = "secs")))

  message(sprintf("  → getStats total: %.2fs", difftime(Sys.time(), t_start, units = "secs")))
  return(stats)
}

#' Extract sparse matrix attributes to mathJS-like sparse matrix format
#'
#' @param matrix sparse matrix
#'
#' @return list with sparse matrix attributes
#'
toSparseJson <- function(matrix) {
  t_start <- Sys.time()
  
  response_unsafe <- list(
    values = matrix@x,
    index = matrix@i,
    ptr = matrix@p,
    size = matrix@Dim
  )

  t_json_start <- Sys.time()
  response <- lapply(response_unsafe, ensure_is_list_in_json)
  message(sprintf("  ✨ toSparseJson: %.2fs (values: %d, indices: %d)", difftime(Sys.time(), t_json_start, units = "secs"), length(matrix@x), length(matrix@i)))

  return(response)
}
