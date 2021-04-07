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
        # Embedding was computed in data processing. Just to avoid possible change in dimensional reduction technique, 
        # I leave this code to check if it has not been computed. 
        if(!"tsne"%in%names(data@reductions))
            data <- RunTSNE(data,
                        reduction = 'pca', 
                        seed.use = 1,
                        dims = 1:pca_nPCs, 
                        perplexity = config$perplexity, 
                        learning.rate = config$learningRate)

        df_embedding <- Embeddings(data, reduction = type)
    } else if(type=="umap"){
        # Embedding was computed in data processing. Just to avoid possible change in dimensional reduction technique, 
        # I leave this code to check if it has not been computed. 
        if(!"umap"%in%names(data@reductions))
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
    
    # Check the size of the cellset to create an empty data.frame. Fill the dataframe with NAs value
    max_cells_id <- length(data$cells_id)-1
    df_embedding_whole <- data.frame(Var_1=rep(NA, max_cells_id), Var_2=rep(NA, max_cells_id))
    # Rename the colnames
    colnames(df_embedding_whole) <- colnames(df_embedding)

    # Index the rows by cells id
    rownames(df_embedding_whole) <- 0:(max_cells_id-1)
    # Get ordered the current cells id of the experiment
    cells_id <- data@meta.data[order(data@meta.data$cells_id), "cells_id"]
    # Update the unfiltered cells id with the embedding information
    df_embedding_whole[as.character(cells_id), ] <- df_embedding

    # Convert to matrix to keep the structure for the API [[x1, y1], [x2, y2] [null, null], ...]
    return(as.matrix(df_embedding_whole))

}