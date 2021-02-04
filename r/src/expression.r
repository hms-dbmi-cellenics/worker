#
# runExpression  
# returns the expression values for each gene in the input.
# 
# req$body has
# genes: list of gene common names to search for in the annotation.  
#
runExpression <- function(req) {
    #Get the annotation matrix with the geneid to name translation, and the subset with the correct names.
    df <- data@misc$gene_annotations
    genesSubset <- subset(df,df$name %in% req$body$genes)
    #Get the expression values for those genes in the corresponding matrix.
    geneExpression <- as.matrix(data@assays$RNA@data[genesSubset$input,])
    #
    #There might be a better way to do this, but when the count is 1 the previous
    #search returns a transposed matrix without the gene id corresponding to that only column
    #so I had to fix that.
    #
    if (ncol(geneExpression)==1){
        geneExpression<-t(geneExpression)
        for(gene in genesSubset$input){
            if (gene %in% rownames(data@assays$RNA@data)){
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
    newNames = c()
    for (gene in rownames(geneExpression)){ 
        newNames<-append(newNames,genesSubset[genesSubset$input == gene,]$name)
    }
    rownames(geneExpression) <- newNames
    geneExpression <- as.list(as.data.frame(t(as.matrix(geneExpression))))
    return(geneExpression)
}

