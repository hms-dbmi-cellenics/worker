#
#runEmbedding 
#Returns a list of x,y coordinates for each cell.
#req is the request. 
#
#Req$body has:
#type = type of embedding, supported umap, pca and tsne.  
#config = config list. 
#
#Config has=
#UMAP:
#minimumDistance = float
#distanceMetric = string (euclidean, cosine, etc)
#
#tsne:
#perplexity
#lerarningRate
#
runEmbedding <- function(req) {
    type <- req$body$type
    config <- req$body$config
    pca_nPCs <- 30 
    if (type == "pca") {
        # Leaving this here to add parameters in the future. Won't leave uncommented to avoid recalculating PCA>
        # RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)
        df_embedding <- Embeddings(data, reduction = type)[,1:2]
    } else if(type=="tsne"){
        data <- RunTSNE(data,
                        reduction = 'pca', 
                        seed.use = 1,
                        dims = 1:pca_nPCs, 
                        perplexity = config$perplexity, 
                        learning.rate = config$learningRate)
        df_embedding <- Embeddings(data, reduction = type)
    } else if(type=="umap"){
        data <- RunUMAP(data,
                        seed.use = 42,
                        reduction='pca', 
                        dims = 1:pca_nPCs, 
                        verbose = F, 
                        min.dist = config$minimumDistance, 
                        metric = config$distanceMetric,
                        umap.method = "uwot-learn")                        
        df_embedding <- Embeddings(data, reduction = type)
    }
    # Order embedding by cells id in ascending form

    df_embedding <- df_embedding[rownames(data@meta.data[order(data@meta.data$cells_id), ]), ]
    ids <- data@meta.data$cells_id+1
    idmax <- max(ids)
    ridx <- seq_len(idmax)
    res <- data.frame(row = ridx, dim1 = NA, dim2 = NA)
    res[ids, ]$dim1 <- df_embedding[,1]
    res[ids, ]$dim2 <- df_embedding[,2]
    res$row <- NULL
    # convert to list of vectors
    na.res <- is.na(res$dim1)
    res <- data.frame(t(res))
    res <- as.list(res)
    res[na.res] <- NA
    res <- unname(res)
    return(res)
}