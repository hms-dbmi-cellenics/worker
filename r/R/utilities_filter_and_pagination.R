applyFilters <- function(gene_results, filters) {
  filter_columns <- sapply(filters, `[[`, "columnName")

  # apply gene text filter
  gene.idx <- which(filter_columns == "gene_names")[1]
  if (!is.na(gene.idx)) {
    gene <- filters[[gene.idx]]$expression
    gene_results <- gene_results[grepl(gene, gene_results$gene_names, ignore.case = TRUE), ]
  }

  # apply numeric filters
  numeric_columns <- colnames(gene_results)[sapply(gene_results, is.numeric)]

  for (idx in seq_along(filter_columns)) {
    column <- filter_columns[idx]
    if (!column %in% numeric_columns) next()
    filter <- filters[[idx]]

    gene_results <- applyNumericFilter(gene_results, column, filter$comparison, filter$value)
  }

  return(gene_results)
}

applyNumericFilter <- function(gene_results, column, comparison, value) {
  if (comparison == "greaterThan") {
    keep <- gene_results[[column]] > value
  } else if (comparison == "lessThan") {
    keep <- gene_results[[column]] < value
  }

  gene_results <- gene_results[keep, ]
  return(gene_results)
}

handlePagination <- function(gene_results, offset, limit, order_by, order_decreasing) {
  full_count <- nrow(gene_results)

  if (order_by %in% names(gene_results)) {
    gene_results <- gene_results[order(gene_results[, order_by], decreasing = order_decreasing), ]
  }

  # Offset starts at 0, limit is number of results per page
  offset <- offset + 1
  limit <- limit - 1

  gene_results <- na.omit(gene_results[(offset):(offset + limit), ])
  return(list(gene_results = gene_results, full_count = full_count))
}
