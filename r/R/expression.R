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
    ignore.case <- TRUE
    if (ignore.case){
        genesSubset <- subset(df, toupper(df$name) %in% toupper(req$body$genes))
    }else{
        genesSubset <- subset(df, df$name %in% req$body$genes)
    }
    
    #
    #Get the expression values for those genes in the corresponding matrix.
    geneExpression <-data@assays$RNA@data[unique(genesSubset$input),,drop=FALSE]
    geneExpression <- as.data.frame(t(as.matrix(geneExpression)))
    geneExpression$cells_id <- data@meta.data$cells_id
    geneExpression <- geneExpression[ order(geneExpression$cells_id), ]
    geneExpression <- geneExpression %>%
        tidyr::complete(cells_id = seq(0,max(data@meta.data$cells_id))) %>%
        select(-cells_id)
    # worried about duplicate gene row.names in @data
    symbol_idx <- match(colnames(geneExpression), genesSubset$input)
    colnames(geneExpression) <- genesSubset$name[symbol_idx]
    #geneExpression <- as.list(as.data.frame(t(as.matrix(geneExpression))))
    return(geneExpression)
}
