#
#
# req has body with        {'selectFields'= {"gene_names", "dispersions"},
#                          'orderBy'= {"gene_names", "dispersions"},
#                          'geneNamesFilter' = string with r placeholders (^ before for endswith and $ after for startswith),
#                          'orderDirection'= {DESC, ASC}
#                          'offset'= int: number of genes from the first one in the sorted and filtered dataframe to offset the results by
#                          'limit'= int: how many results per page.
#
# returns: df with variance.standarized, SYMBOL (gene_names) and full_count collumn with number of genes after filtering.
#
# We decided to return variance.standarized because it's a commonly used indicator of dispersion that accounts for the mean.[1]
#
# [1]Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.
#
#' @export
getList <- function(req, data) {
  select_fields <- req$body$selectFields
  order_by <- req$body$orderBy
  order_decreasing <- req$body$orderDirection == "DESC"
  offset <- req$body$offset
  limit <- req$body$limit

  # Gene dispersion slot generated in data ingest script with the same info as the meta.features slot but with the annotated genes
  gene_results <- data@misc$gene_dispersion
  colnames(gene_results)[1:3] <- c('mean', 'variance', 'variance.standardized')

  colnames(gene_results)[colnames(gene_results) == "SYMBOL"] <- "gene_names"
  colnames(gene_results)[colnames(gene_results) == "variance.standardized"] <- "dispersions"

  # apply gene name filter
  gene_pattern <- req$body$geneNamesFilter
  if (!is.null(gene_pattern)) {
    gene_filter <- list(list(columnName = "gene_names", expression = gene_pattern))
    gene_results <- applyFilters(gene_results, gene_filter)
  }

  paginated_results <- handlePagination(gene_results, offset, limit, order_by, order_decreasing)
  gene_results <- paginated_results$gene_results

  columns <- c("gene_names", "dispersions")

  gene_results <- gene_results[, columns]

  return(list(gene_results = gene_results, full_count = paginated_results$full_count))
}
