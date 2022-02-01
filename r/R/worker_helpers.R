
#
# subset_ids subsets a seurat object with the cell ids
#
subset_ids <- function(scdata, cells_id) {
  meta_data_subset <- scdata@meta.data[match(cells_id, scdata@meta.data$cells_id), ]
  current_cells <- rownames(meta_data_subset)
  scdata <- subset(scdata, cells = current_cells)
  return(scdata)
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
    dplyr::slice(which.min(pval))

  top_markers <- filtered_markers %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::arrange(group)

  message(sprintf("%d markers selected", nrow(top_markers)))
  return(top_markers)
}

getMarkerNames <- function(data, all_markers) {
  all_markers$name <- data@misc$gene_annotations[all_markers$feature, "name"]
  all_markers <- all_markers %>% dplyr::transmute(input = feature, name = name)
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
    dplyr::select(-cells_id)
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

applyFilters <- function(gene_results, filters) {
  filter_columns <- sapply(filters, `[[`, "columnName")

  # apply gene text filter
  gene.idx <- which(filter_columns == "gene_names")[1]
  if (!is.na(gene.idx)) {
    gene <- filters[[gene.idx]]$expression
    gene_results <- gene_results[grepl(gene, gene_results$gene_names, ignore.case = TRUE), ]
  }

  # apply numeric filters
  numeric_columns <- c("logFC", "p_val_adj", "pct_1", "pct_2", "auc")

  for (idx in seq_along(filter_columns)) {
    column <- filter_columns[idx]
    if (!column %in% numeric_columns) next()
    filter <- filters[[idx]]

    gene_results <- applyNumericFilter(gene_results, column, filter$comparison, filter$value)
  }

  return(gene_results)
}

applyNumericFilter <- function(gene_results, column, comparison, value) {
  if (comparison == "greaterThan") {
    keep <- gene_results[[column]] > value
  } else if (comparison == "lessThan") {
    keep <- gene_results[[column]] < value
  }

  gene_results <- gene_results[keep, ]
  return(gene_results)
}

handlePagination <- function(gene_results, offset, limit, order_by, order_decreasing) {
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

insert_set_child_through_api <- function(new_cell_set, api_url, experiment_id, cell_set_key, auth_JWT) {
  httr_query <- paste0("$[?(@.key == \"", cell_set_key, "\")]")
  children <- list("$insert" = list(index = "-", value = new_cell_set))

  httr::PATCH(
    paste0(api_url, "/v1/experiments/", experiment_id, "/cellSets"),
    body = list(
      list("$match" = list(
        query = httr_query,
        value = list("children" = children)
      )),
      encode = "json",
      httr::add_headers(
        "Content-Type" = "application/boschni-json-merger+json",
        "Authorization" = auth_JWT
      )
    )
  )
}
