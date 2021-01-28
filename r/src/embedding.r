runEmbedding <- function (req,type) {
    if(type=="pca"){
        #Leaving this here to add parameters in the future.
        #RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)
        return(Embeddings(data, reduction = "pca"))
    }else{
        #data <- FindNeighbors(data, k.param = 20, annoy.metric = "cosine", verbose=FALSE) #default method
        #message("Running embedding")
        #data <- RunUMAP(data, reduction='pca', dims = 1:10, verbose = F, umap.method = "uwot-learn")
        return(Embeddings(data, reduction = "umap"))
    }
}
