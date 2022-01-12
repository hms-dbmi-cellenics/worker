#
# runExpression
# returns the expression values for each gene in the input.
#
# req$body has
# genes: list of gene common names to search for in the annotation.
#
#
#For now we return the values stored in data (normalized values). When the correct config parameter is set on the UI, we'll add the scaled values.
#
#' @export
runExpression <- function(req, data) {
    #Get the annotation matrix with the geneid to name translation, and the subset with the correct names.
    df <- data@misc$gene_annotations
    genesSubset <- subset(df, toupper(df$name) %in% toupper(req$body$genes))
    genesSubset <- genesSubset[,c("input","name")]
    return(list())
}
