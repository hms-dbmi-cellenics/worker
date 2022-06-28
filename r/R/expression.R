#
# runExpression
# returns the expression values for each gene in the input.
#
# req$body has
# genes: list of gene common names to search for in the annotation.
#
#
# For now we return the values stored in data (normalized values). When the correct config parameter is set on the UI, we'll add the scaled values.
#
#' @export
runExpression <- function(req, data) {
  # Get the annotation matrix with the geneid to name translation, and the subset with the correct names.
  df <- data@misc$gene_annotations
  genesSubset <- subset(df, toupper(df$name) %in% toupper(req$body$genes))

  if (!nrow(genesSubset)) {
    stop(
      generateErrorMessage(
        error_codes$GENE_NOT_FOUND,
        paste("Gene(s):", paste(req$body$genes, collapse = ", "), "not found!")
      )
    )
  }

  genesSubset <- genesSubset[, c("input", "name")]
  res <- getExpressionValues(genesSubset, data)
  res <- formatExpression(res)
  return(res)
}
