
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
