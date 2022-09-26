# runExpression
# returns the expression values for each gene in the input.
#
# req$body has
# genes: list of gene common names to search for in the annotation.
#
#' @export
runExpression <- function(req, data) {
  gene_annotations <- data@misc$gene_annotations

  # subset with gene NAMES passed from UI
  gene_subset <-
    subset(
      gene_annotations,
      toupper(gene_annotations$name) %in% toupper(req$body$genes)
    )

  if (!nrow(gene_subset)) {
    stop(generateErrorMessage(
      error_codes$GENE_NOT_FOUND,
      paste(
        "Gene(s):",
        paste(req$body$genes, collapse = ", "),
        "not found!"
      )
    ))
  }

  gene_subset <- gene_subset[, c("input", "name")]

  return(getGeneExpression(data, gene_subset))
}
