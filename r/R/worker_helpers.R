
#
# subset_ids subsets a seurat object with the cell ids
#
subset_ids <- function(scdata, cells_id) {
  meta_data_subset <- scdata@meta.data[match(cells_id, scdata@meta.data$cells_id), ]
  current_cells <- rownames(meta_data_subset)
  scdata <- subset(scdata, cells = current_cells)
  return(scdata)
}

get_incomplete_clusters <- function(top_markers, nFeatures) {
  top_markers %>%
    dplyr::group_by(group) %>%
    dplyr::tally() %>%
    dplyr::filter(n < nFeatures) %>%
    dplyr::mutate(missing_genes = nFeatures - n)
}

unique_missing_markers <- function(all_markers, top_markers, clusters) {
  # returns tibble with markers of incomplete clusters but not present in
  # top markers
  all_markers %>%
    # remove markers present in other clusters
    dplyr::filter(group %in% clusters$group) %>%
    dplyr::anti_join(top_markers, by = "feature") %>%
    # get number of genes needed
    dplyr::left_join(clusters, by = "group") %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(.data$pval)) %>%
    dplyr::ungroup()
}

get_low_quality_markers <- function(all_markers, top_markers, clusters, nFeatures) {
  # returns tibble with new markers of incomplete clusters

  unique_markers <- unique_missing_markers(all_markers, top_markers, clusters)
  message(sprintf("unique_markers: %d", nrow(unique_markers)))

  unique_markers %>%
    # use logFC as decision criteria
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::group_by(group) %>%
    # actually get markers
    dplyr::slice(seq(dplyr::first(missing_genes))) %>%
    dplyr::select(-missing_genes) %>%
    dplyr::ungroup()
}

getTopMarkerGenes <- function(nFeatures, data, cellSets, aucMin = 0.3, pctInMin = 20, pctOutMax = 70) {
  data$marker_groups <- NA

  object_ids <- data$cells_id
  for (i in seq_along(cellSets)) {
    set <- cellSets[[i]]
    filtered_cells <- intersect(set$cellIds, object_ids)
    data$marker_groups[object_ids %in% filtered_cells] <- i
  }

  all_markers <- presto::wilcoxauc(data, group_by = "marker_groups", assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)

  # may not return nFeatures markers per cluster if values are too stringent
  filtered_markers <- all_markers %>%
    dplyr::filter(logFC > 0 &
      auc >= aucMin &
      pct_in >= pctInMin &
      pct_out <= pctOutMax) %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(pval)) %>%
    ungroup()

  top_markers <- filtered_markers %>%
    dplyr::arrange(group, desc(logFC)) %>%
    dplyr::group_by(group) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::ungroup()

  # check if there are incomplete clusters
  incomplete_clusters <- get_incomplete_clusters(top_markers, nFeatures)
  message("Incomplete clusters:")
  message(print(incomplete_clusters))

  if (nrow(incomplete_clusters > 0)) {
    extra_markers <- get_low_quality_markers(
      all_markers,
      top_markers,
      incomplete_clusters,
      nFeatures
    )
    top_markers <- dplyr::bind_rows(top_markers, extra_markers) %>% dplyr::arrange(group, dplyr::desc(logFC))
  }

  message(sprintf("%d markers selected", nrow(top_markers)))
  return(top_markers)
}

getMarkerNames <- function(data, all_markers) {
  all_markers$name <- data@misc$gene_annotations[all_markers$feature, "name"]
  all_markers <- all_markers %>% transmute(input = feature, name = name)
  rownames(all_markers) <- c()
  return(all_markers)
}

#' Returns expression values for selected genes
#'
#' @param genes - Must have names and input(ensmbl ids)
#' @param data
#'
#' @return
#' @export
#'
#' @examples
getExpressionValues <- function(genes, data) {
  quantile_threshold <- 0.95
  #
  # Get the expression values for those genes in the corresponding matrix.
  geneExpression <- data@assays$RNA@data[unique(genes$input), , drop = FALSE]
  geneExpression <- as.data.frame(t(as.matrix(geneExpression)))
  geneExpression$cells_id <- data@meta.data$cells_id
  geneExpression <- geneExpression[order(geneExpression$cells_id), ]
  geneExpression <- geneExpression %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    select(-cells_id)
  # worried about duplicate gene row.names in @data
  symbol_idx <- match(colnames(geneExpression), genes$input)
  colnames(geneExpression) <- genes$name[symbol_idx]
  adjGeneExpression <- as.data.frame(apply(geneExpression, 2, FUN = function(x) {
    lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
    i <- 0.01
    while (lim == 0 & i + quantile_threshold <= 1) {
      lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
      i <- i + 0.01
    }
    return(pmin(x, lim))
  }))
  return(list(rawExpression = geneExpression, truncatedExpression = adjGeneExpression))
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
  snn_name <- paste0(data@active.assay, "_snn")

  # if doesn't exist, run SNN
  if (!snn_name %in% names(data)) data <- Seurat::FindNeighbors(data)

  # convert Seurat Graph object to igraph
  # similar to https://github.com/joshpeters/westerlund/blob/46609a68855d64ed06f436a6e2628578248d3237/R/functions.R#L85
  adj_matrix <- Matrix::Matrix(as.matrix(data@graphs[[snn_name]]), sparse = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE
  )


  return(g)
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
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resolution))
  algorithm <- list("louvain" = 1, "leiden" = 4)[[type]]
  # To run clustering, we need to identify the active.reduction that is used in FindNeighbors.
  if ("active.reduction" %in% names(data@misc)) {
    active.reduction <- data@misc[["active.reduction"]]
  } else {
    active.reduction <- "pca"
  }

  if (type == "leiden") {
    # emulate FindClusters, which overwrites seurat_clusters slot and meta.data column
    g <- getSNNiGraph(data)
    clus_res <- igraph::cluster_leiden(g, "modularity", resolution_parameter = resolution)
    clusters <- clus_res$membership
    names(clusters) <- clus_res$names
    clusters <- clusters[colnames(data)]
    data$seurat_clusters <- data@meta.data[, res_col] <- factor(clusters - 1)
  } else {
    graph.name <- paste0(DefaultAssay(data), "_snn")
    if (!graph.name %in% names(data)) {
      data <- Seurat::FindNeighbors(data, k.param = 20, annoy.metric = "cosine", verbose = FALSE, reduction = active.reduction)
    }
    data <- FindClusters(data, resolution = resolution, verbose = FALSE, algorithm = algorithm)
  }
  return(data)
}

handle_pagination <- function(gene_results, offset, limit, order_by, order_decreasing, filter) {
  if (!is.null(filter)) {
    gene_results <- gene_results[grepl(filter, gene_results$gene_name, ignore.case = TRUE), ]
  }

  full_count <- nrow(gene_results)

  if (order_by %in% names(gene_results)) {
    gene_results <- gene_results[order(gene_results[, order_by], decreasing = order_decreasing), ]
  }

  # Offset starts at 0, limit is number of results per page
  offset <- offset + 1
  limit <- limit - 1

  gene_results <- na.omit(gene_results[(offset):(offset + limit), ])
  return(list(gene_results = gene_results, full_count = full_count))
}
