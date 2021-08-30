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
    df <- data@misc$gene_annotations
    row.names(df) <- df$name

    genes <- req$body$genes
    enids <- df[genes, 'input']

    expr <- data[['RNA']]@data[enids,, drop = FALSE] %>% as.matrix()

    raw_expr <- expr %>%
        completeCellIds(data$cells_id) %>%
        `colnames<-`(genes)

    trunc_expr <- expr %>%
        truncateExpression() %>%
        completeCellIds(data$cells_id) %>%
        `colnames<-`(genes)

    res <- list(
        rawExpression = raw_expr,
        truncatedExpression = trunc_expr
    )
    return(res)
}

runHeatmapExpression <- function(req, data) {
    df <- data@misc$gene_annotations
    row.names(df) <- df$name

    enids <- df[req$body$genes, 'input']

    getHeatmapExpression(data, enids)
}
