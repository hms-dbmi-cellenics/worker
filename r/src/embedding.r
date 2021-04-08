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
    ids <- data@meta.data$cells_id
    idmax <- max(ids)
    ids <- ids+1
    tmp <- data.frame(row = 0:idmax)
    tmp$UMAP_1 <- tmp$UMAP_2 <- NA
    tmp[ids, ]$UMAP_1 <- df_embedding[,1]
    tmp[ids, ]$UMAP_2 <- df_embedding[,2]
    tmp$row <- NULL
    row.names(tmp)[ids] <- row.names(df_embedding)
    tmp <- lapply(row.names(tmp), function(rname) {
        umap1 <- tmp[rname, 'UMAP_1']
        umap2 <- tmp[rname, 'UMAP_2']
         if(is.na(umap1)) {
            return(NA) 
         } else {
            return(c(umap1, umap2))
         }
    })
    return(tmp)


    # Order embedding by cells id in ascending form
    #df_embedding <- df_embedding[rownames(data@meta.data[order(data@meta.data$cells_id), ]), ]

    #return(df_embedding)

}