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
  type <- req$body$type
  algo <- list("louvain"=1,"leiden"=4)[[type]]
  res_col <- paste0(data@active.assay, "_snn_res.",toString(resol))
  
  if (type == 'leiden') {
    # emulate FindClusters, which overwrites seurat_clusters slot and meta.data column
    g <- getSNNiGraph(data)
    clus_res <- igraph::cluster_leiden(g, 'modularity', resolution_parameter = resol)
    clusters <- clus_res$membership
    names(clusters) <- clus_res$names
    clusters <- clusters[colnames(data)]
    data$seurat_clusters <- data@meta.data[, res_col] <- factor(clusters-1)
    
  } else {
    data <- FindClusters(data, resolution=resol, verbose = FALSE, algorithm = algo) 
  }
  
  #In the meta data slot the clustering is stored with the resolution used to calculate it
  # RNA_snn_res.#resolution
  df <- data.frame(cluster=data@meta.data[,res_col], cell_ids=data@meta.data$cells_id)  
  #get the cell barcodes as rownames
  rownames(df) <- rownames(data@meta.data)
  return(df)
}


#' Get and Convert SNN Graph object into igraph object
#' 
#' This is used to facilitate leiden clustering.
#'
#' @param data \code{Seurat} object
#'
#' @return boolean indicating if SNN Graph object exists
#' 
getSNNiGraph <- function(data) {
  
  # check to see if we already have Seurat SNN Graph object
  snn_name <- paste0(data@active.assay, '_snn')
  
  # if doesn't exist, run SNN
  if (!snn_name %in% names(data)) data <- Seurat::FindNeighbors(data)
  
  # convert Seurat Graph object to igraph
  # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
  adj_matrix <- Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, 
                                           mode = 'undirected',
                                           weighted = TRUE)
  
  
  return(g)
}

