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
runExpression <- function(req) {
    #Get the annotation matrix with the geneid to name translation, and the subset with the correct names.
    df <- data@misc$gene_annotations
    genesSubset <- subset(df,df$name %in% req$body$genes)
    #Get the expression values for those genes in the corresponding matrix.
    geneExpression <- as.matrix(data@assays[[data@active.assay]]@data[genesSubset$input,])
    #
    #There might be a better way to do this, but when the count is 1 the previous
    #search returns a transposed matrix without the gene id corresponding to that only column
    #so I had to fix that.
    #
    if (ncol(geneExpression)==1){
        geneExpression<-t(geneExpression)
        for(gene in genesSubset$input){
            #
            #data@assays$RNA@data
            if (gene %in% rownames(data@assays[[data@active.assay]]@data)){
                rownames(geneExpression)<-gene  
                break
            }    
        }
    }
    #
    #Now I need to get the gene name from the IDs of the ones that were found on the data matrix.
    #With those values, replace the rownames. If I had used the whole genesSubset$input as gene names
    #I'd have risked adding a column for a gene that wasn't found on the data matrix.
    #
    #Another possible way of doing this is:
    #merge(genesSubset,geneExpression,by.x="input",by.y=0,all.y=TRUE,all.x=FALSE)
    #I don't know which one is faster, for now I'll leave this one.
    #
    #Juanlu suggested genesSubset$name[match(rownames(geneExpression), genesSubset$input)]
    #
    rownames(geneExpression) <- genesSubset$name[match(rownames(geneExpression), genesSubset$input)]
    geneExpression <- as.list(as.data.frame(t(as.matrix(geneExpression))))
    return(geneExpression)
}

