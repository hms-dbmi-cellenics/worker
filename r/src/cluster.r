#
# getClusters 
# returns the clusters in the shape of a dataframe with a clusters column, cell ids column and cell barcode as rownames. 
#
# req$body has:
# type: can be "louvain"/"leiden"
# config:{
#          resolution: integer, range: 0 - 2         
#         }
#
getClusters <- function(req){
    resol <- req$body$config$resolution
    algo <- list("louvain"=1,"leiden"=4)[[req$body$type]]
    #Leaving neighbors here in case we eventually set parameters.
    #data <- FindNeighbors(data, k.param = 20, annoy.metric = neighbors_metric, verbose=FALSE)
    data <- FindClusters(data, resolution=resol, verbose = FALSE, algorithm = algo) 
    #In the meta data slot the clustering is stored with the resolution used to calculate it
    # RNA_snn_res.#resolution
    #

    str <- paste(data@active.assay,"_snn_res.",toString(resol),sep = "")
    df <- data.frame("cluster"= data@meta.data[,str], "cell_ids"=data@meta.data$cells_id)  
    #get the cell barcodes as rownames
    rownames(df) <- rownames(data@meta.data)
    return(df)
}
