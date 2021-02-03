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
runEmbedding <- function(req) {
    type <- req$body$type
    config <- req$body$config
    if (type == "pca") {
        # Leaving this here to add parameters in the future. Won't leave uncommented to avoid recalculating PCA>
        # RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)
        return(Embeddings(data, reduction = type)[,1:2])
    } else if(type=="tsne"){
        # Need to add parameters to the UI for tsne.
        data <- RunTSNE(data,
                        reduction="pca",
                        dims = 1:10                      
                        ) 
        return(Embeddings(data, reduction = type))
    } else if(type=="umap"){
        #Checking that the parameters are correctly set. 
        if ("minimumDistance" %in% names(config)){
            mindist <- config$minimumDistance
        }else{
            mindist <- 0.3
        }
        if ("distanceMetric" %in% names(config)){
            umapMetric <- config$distanceMetric
        }else{
            umapMetric <- "euclidean"
        }
        data <- RunUMAP(data,
                        reduction='pca', 
                        dims = 1:10, 
                        verbose = F, 
                        min.dist = config$minimumDistance, 
                        metric = config$distanceMetric,
                        umap.method = "uwot-learn")                        
        return(Embeddings(data, reduction = type))        
    }
}