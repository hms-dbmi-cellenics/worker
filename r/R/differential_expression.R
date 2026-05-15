#' Run Differential Expression Task
#'
#' Runs either a comparison within samples/groups (i.e. marker genes) or
#' between samples/groups (differential expression analysis). The request post
#' includes the cell IDs that we need to compare in a cell set named "base"
#' and "background".
#'
#' @export
#'
runDE <- function(req, data) {

  # add comparison group to 'custom' slot
  data <- addComparisonGroup(req, data)

  comparison_type <- req$body$comparisonType
  if (comparison_type == "within") {
    result <- runWilcoxAUC(data)
  } else if (comparison_type == "between") {
    pbulk <- makePseudobulkMatrix(data)
    result <- runPseudobulkDE(pbulk)
  }

  # replace name with gene names and add Ensembl IDs
  annot <- data@misc$gene_annotations[row.names(result),]
  result$gene_names <- annot$name
  result$Gene <- annot$input

  # replace NA gene symbols with ensembl ids
  na_genes <- is.na(result$gene_names)
  result$gene_names[na_genes] <- result$Gene[na_genes]

  # replace 0 in p_val_adj with the smallest floating-point value
  # this is required to correctly plot log(p_val_adj)
  # in the volcano plot, because log(0)=Inf
  if ("p_val_adj" %in% names(result)) {
    result["p_val_adj"][result["p_val_adj"] == 0] <- .Machine$double.xmin
  }

  if ("pagination" %in% names(req$body)) {
    result <- paginateDE(result, req)
  } else {
    result <- list(gene_results = result, full_count = nrow(result))
  }

  return(result)
}

runWilcoxAUC <- function(data) {

  # for speed: take at most 10000 cells per group
  set.seed(0)
  max_cells <- 10000

  X_matrix <- data[["RNA"]]$data
  y <- data$custom

  keep_cell_ids <- data@meta.data |>
    dplyr::group_by(custom) |>
    dplyr::slice_sample(n = max_cells, replace = FALSE) |>
    dplyr::ungroup() |>
    dplyr::pull(cells_id)

  if (length(keep_cell_ids) < length(y))
    message(
      "Using at most 10,000 cells per group for marker gene testing.",
      " Total cells used: ", length(keep_cell_ids), "."
    )

  object_ids <- data$cells_id
  keep_indices <- match(keep_cell_ids, object_ids)
  mat_subset <- X_matrix[, keep_indices]
  mat_subset <- as(mat_subset, "dgCMatrix")
  group_subset <- y[keep_indices]

  result <- presto::wilcoxauc(mat_subset, y = group_subset)
  result <- result[result$group == "base", ]

  rownames(result) <- result$feature
  result <- result[, c("pval", "logFC", "pct_in", "pct_out", "padj", "auc")]
  colnames(result) <- c("p_val", "logFC", "pct_1", "pct_2", "p_val_adj", "auc")

  return(result)
}

paginateDE <- function(result, req) {
  message("Paginating results:  ", str(result))
  pagination <- req$body$pagination
  genes_only <- FALSE

  order_by <- pagination$orderBy
  order_decreasing <- pagination$orderDirection == "DESC"
  offset <- pagination$offset
  limit <- pagination$limit
  filters <- pagination$filters

  result <- applyFilters(result, filters)

  if ("genesOnly" %in% names(req$body)) {
    genes_only <- req$body$genesOnly
  }

  if (genes_only) {
    n_genes <- min(limit, nrow(result))

    result <- result[order(result[, order_by], decreasing = order_decreasing), ]

    gene_results <- list(
      gene_names = result$gene_names[1:n_genes],
      gene_id = result$Gene[1:n_genes]
    )

    result <- list(gene_results = gene_results, full_count = n_genes)
    return(result)
  }

  result <- handlePagination(result, offset, limit, order_by, order_decreasing)

  return(result)
}
