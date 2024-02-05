# IMPORTANT: functions in this file are duplicated in the pipeline.
# If you update, change both.

#' Get Clusters
#'
#' @param data SeuratObject
#' @param req list with items:
#'  \code{req$body$type} either "louvain" or "leiden"
#'  \code{req$body$config$resolution} integer, range: 0 - 2
#'
#' @return data.frame with columns "cluster" and "cell_ids"
#' @export
#'
#' @examples
runClusters <- function(req, data) {
  type <- req$body$type
  config <- req$body$config
  resolution <- config$resolution


  data <- getClusters(type, resolution, data)
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
  # In the meta data slot the clustering is stored with the resolution
  # used to calculate it
  # RNA_snn_res.#resolution
  df <- data.frame(
    cluster = data@meta.data[, res_col],
    cell_ids = data@meta.data$cells_id
  )
  # get the cell barcodes as rownames
  rownames(df) <- rownames(data@meta.data)

  message("formatting cellsets")
  formatted_cell_sets <-
    format_cell_sets_object(df, type, data@misc$color_pool)

  message("updating through api")
  updateCellSetsThroughApi(
    formatted_cell_sets,
    req$body$apiUrl,
    req$body$experimentId,
    "louvain",
    req$body$authJwt,
    append = FALSE
  )

  return(df)
}

#' Compute clusters and return object with clusters
#'
#' @param algorithm
#' @param resolution
#' @param data
#'
#' @return
#'
#' @examples
getClusters <- function(type, resolution, data) {
  res_col <-
    paste0(data@active.assay, "_snn_res.", toString(resolution))
  algorithm <- list("louvain" = 1, "leiden" = 4)[[type]]

  # use the reduction from data integration for nearest neighbors graph
  if ("active.reduction" %in% names(data@misc)) {
    active.reduction <- data@misc[["active.reduction"]]
  } else {
    active.reduction <- "pca"
  }

  if (type == "leiden") {
    # emulate FindClusters, which overwrites seurat_clusters slot and meta.data
    # column
    snn_graph <- getSNNiGraph(data, active.reduction)
    clus_res <-
      igraph::cluster_leiden(snn_graph, "modularity", resolution_parameter = resolution)
    clusters <- clus_res$membership
    names(clusters) <- clus_res$names
    clusters <- clusters[colnames(data)]
    data$seurat_clusters <-
      data@meta.data[, res_col] <- factor(clusters - 1)
  } else {
    graph_name <- paste0(Seurat::DefaultAssay(data), "_snn")
    if (!graph_name %in% names(data)) {
      # number of dimensions used must be lte to available dimensions
      dims <- 1:min(10, length(data@reductions[[active.reduction]]))
      data <-
        Seurat::FindNeighbors(
          data,
          annoy.metric = "cosine",
          reduction = active.reduction,
          dims = dims,
          verbose = FALSE,
        )
    }
    data <- Seurat::FindClusters(
      data,
      resolution = resolution,
      verbose = FALSE,
      algorithm = algorithm,
      random.seed = ULTIMATE_SEED
    )
  }

  return(data)
}


#' Get and Convert SNN Graph object into igraph object
#'
#' This is used to facilitate leiden clustering.
#'
#' @param data \code{Seurat} object
#'
#' @return boolean indicating if SNN Graph object exists
#'
getSNNiGraph <- function(data, active.reduction) {
  # check to see if we already have Seurat SNN Graph object
  snn_name <- paste0(data@active.assay, "_snn")

  # if doesn't exist, run SNN
  if (!snn_name %in% names(data)) {
    dims <- 1:min(10, length(data@reductions[[active.reduction]]))
    data <-
      Seurat::FindNeighbors(data,
                            reduction = active.reduction,
                            dims = dims,
                            verbose = FALSE)
  }

  # convert Seurat Graph object to igraph
  # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
  adj_matrix <-
    Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE
  )
  return(graph)
}

#' Formats cell sets object for patching through the API
#'
#' This function is only used to format clustering cellsets. Converting from
#' data.frame to list and adding slots necessary for the cellsets file.
#'
#' @param cell_sets data.frame with two columns: cluster and cell_ids
#' @param clustering_method string Either louvain or leiden.
#' @param color_pool character vector of colors in hex
#'
#' @return list
#' @export
#'
#' @examples
format_cell_sets_object <-
  function(cell_sets, clustering_method, color_pool) {
    name <- paste0(clustering_method, " clusters")
    cell_sets_object <-
      list(
        key = "louvain",
        name = name,
        rootNode = TRUE,
        type = "cellSets",
        children = list()
      )
    for (i in sort(unique(cell_sets$cluster))) {
      cells <- cell_sets[cell_sets$cluster == i, "cell_ids"]
      new_set <- list(
        key = paste0(clustering_method, "-", i),
        name = paste0("Cluster ", i),
        rootNode = FALSE,
        type = "cellSets",
        color = color_pool[1],
        cellIds = unname(cells)
      )
      color_pool <- color_pool[-1]
      cell_sets_object$children <-
        append(cell_sets_object$children, list(new_set))
    }
    return(cell_sets_object)
  }
