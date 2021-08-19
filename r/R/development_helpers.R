#' Generates request for worker task
#'
#' @param task
#'
#' @return
#' @export
#'
#' @examples
generate_request <- function(task) {
    res = switch(task,
           "expression" = list(body=list(genes=list("gzma"))),
           "markerHeatmap" = list(body=list(nGenes=5,type="louvain",config=list(resolution=0.5))),
           "DE" = list(body=list(baseCells=0:500,backgroundCells=501:1499)),
           "umap"=list(body=list(type="umap",config=list(minimumDistance=0.3,distanceMetric="cosine"))),
           "tsne"=list(body=list(type="tsne",config=list(perplexity=1,learningRate=1))),
           "listGenes"=list(body=list(name="ListGenes",selectFields=list("gene_names","dispersions"),orderBy="dispersions",orderDirection="DESC",offset=0,limit=50))
           )
}

#' Loads Experiment data
#'
#' @param exp_id
#'
#' @return
#' @export
#'
#' @examples
load_experiment_data <- function(exp_id="test"){
    return(readRDS(paste0("../data/",exp_id,"/r.rds",sep="")))
}

generate_function_speed_curve <- function(funcion,functionString,rFilesPath="../data/time curve"){
    library(ggplot2)
    req <- generate_request(functionString)
    fileList <- list.files(rFilesPath,full.names=TRUE)
    times <- c()
    experiments <- c()
    files<-c()
    ncells <- c()
    ngenes<-c()
    for(i in fileList){
        r <- readRDS(i)
        start_time <- Sys.time()
        funcion(req,r)
        end_time <- Sys.time()
        times<-append(times,as.numeric(difftime(end_time, start_time, units='secs')))
        experiments<-append(experiments,r@misc$experimentId)
        files<-append(files,i)
        ncells<-append(ncells,ncol(r@assays$RNA@counts))
        ngenes<-append(ngenes,nrow(r@assays$RNA@counts))
    }
    res_df <- data.frame(times)
    res_df$experiments <- experiments
    res_df$files <- files
    res_df$ngenes<- ngenes
    res_df$ncells<-ncells
    rownames(res_df) <- files
    return(res_df)
}
