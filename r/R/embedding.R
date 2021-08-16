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
#' @export
runEmbedding <- function(req, data) {
    type <- req$body$type
    config <- req$body$config
    pca_nPCs <- 30

    # To run embedding, we need to set the reduction.
    if("active.reduction" %in% names(data@misc))
        active.reduction <- data@misc[["active.reduction"]]
    else
        active.reduction <- "pca"

    # The slot numPCs is set in dataIntegration with the selectd PCA by the user.
    if("numPCs" %in% names(data@misc))
        pca_nPCs <- data@misc[["numPCs"]]

    message("Active reduction --> ", active.reduction)
    message("Active numPCs --> ", pca_nPCs)
    message("Number of cells/sample:")
    table(data$samples)

    if (type == "pca") {
        # Leaving this here to add parameters in the future. Won't leave uncommented to avoid recalculating PCA>
        # RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)
        df_embedding <- Embeddings(data, reduction = type)[,1:2]
    } else if(type=="tsne"){
        data <- RunTSNE(data,
                        reduction = active.reduction,
                        seed.use = 1,
                        dims = 1:pca_nPCs,
                        perplexity = config$perplexity,
                        learning.rate = config$learningRate)
        df_embedding <- Embeddings(data, reduction = type)

    } else if (type=="umap") {
        # until we figure out why umap-learn and uwot break locally
        env <- Sys.getenv('CLUSTER_ENV', 'development')
        umap.method <- ifelse(env == 'development', 'uwot-learn', 'umap-learn')
        message(sprintf('CLUSTER_ENV: %s --> UMAP method is: %s', env, umap.method))

        data <- RunUMAP(data,
                        seed.use = 42,
                        reduction=active.reduction,
                        dims = 1:pca_nPCs,
                        verbose = F,
                        min.dist = config$minimumDistance,
                        metric = config$distanceMetric,
                        umap.method = umap.method)

        df_embedding <- Embeddings(data, reduction = type)
    }
    # Order embedding by cells id in ascending form
    df_embedding <- as.data.frame(df_embedding)
    df_embedding$cells_id <- data@meta.data$cells_id
    df_embedding <- df_embedding[ order(df_embedding$cells_id), ]
    df_embedding <- df_embedding %>% tidyr::complete(cells_id = seq(0,max(data@meta.data$cells_id))) %>% select(-cells_id)
    res <- purrr::map2(df_embedding[[1]], df_embedding[[2]], function(x, y) { if(is.na(x)) { return(NULL) } else { return(c(x, y)) } } )
    return(res)
}
