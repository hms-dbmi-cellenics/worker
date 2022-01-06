mock_gene_results <- function() {
  data(exprs, package = "presto")
  data(y, package = "presto")

  result <- presto::wilcoxauc(exprs, y)
  result <- result[result$group == "A", ]
  rownames(result) <- result$feature
  result$gene_names <- result$feature
  result <- result[, c("gene_names", "pval", "logFC", "pct_in", "pct_out", "padj", "auc")]
  colnames(result) <- list(
    "gene_names", "p_val", "logFC", "pct_1", "pct_2", "p_val_adj", "auc"
  )
  return(result)
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )
  return(pbmc_small)
}

test_that("applyFilters works with gene with exact match", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "G5"))

  gene_results <- applyFilters(gene_results, filters)

  expect_equal(gene_results$gene_names, "G5")
})

test_that("applyFilters works with gene with case insensitive match", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "g5"))

  gene_results <- applyFilters(gene_results, filters)

  expect_equal(gene_results$gene_names, "G5")
})

test_that("applyFilters returns all results that match gene search pattern", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "G"))

  filt_results <- applyFilters(gene_results, filters)

  expect_equal(filt_results, gene_results)
})


test_that("applyFilters works with single numeric filter", {
  gene_results <- mock_gene_results()
  filter_up <- list(list(
    columnName = "logFC",
    comparison = "greaterThan", value = 0
  ))
  filter_dn <- list(list(
    columnName = "logFC",
    comparison = "lessThan", value = 0
  ))

  results_up <- applyFilters(gene_results, filter_up)
  results_dn <- applyFilters(gene_results, filter_dn)

  expect_equal(results_up, gene_results[gene_results$logFC > 0, ])
  expect_equal(results_dn, gene_results[gene_results$logFC < 0, ])
})

test_that("running applyFilters sequentially is equivalent to running together", {
  gene_results <- mock_gene_results()
  filter_gene <- list(columnName = "gene_names", expression = "G1")
  filter_up <- list(columnName = "logFC", comparison = "greaterThan", value = 0)

  # apply sequentially
  results_seq <- applyFilters(gene_results, list(filter_gene))
  results_seq <- applyFilters(results_seq, list(filter_up))

  # apply together
  results_together <- applyFilters(gene_results, list(filter_gene, filter_up))

  expect_equal(results_seq, results_together)
})

# subset_ids

test_that("subset_ids subsets seurat object correctly", {})

test_that("subset_ids returns same seurat object if empty cells_id", {})

# getTopMarkerGenes

test_that("getTopMarkerGenes returns an object with correct columns", {})

test_that("getTopMarkerGenes returns n_genes * n_clusters at the most", {})

test_that("getTopMarkerGenes returns correctly filtered marker genes", {})

test_that("getTopMarkerGenes returns empty if no genes match filters", {})

# getMarkerNames

test_that("getMarkerNames returns object with correct columns", {})

test_that("getMarkerNames returns same number of gene names as gene_ids", {})

test_that("getMarkerNames returns same number of gene names as gene_ids", {})

# this might not make sense, the name column contains the ID if there's no gene
# name
test_that("getMarkerNames returns gene_id if there's no gene name", {})


# getExpressionValues

test_that("getExpressionValues returns object with correct elements", {})

test_that("getExpressionValues returns same number of values as genes requested", {})

test_that("getExpressionValues adjusted values are correct", {})

# getSNNigraph

test_that("getSNNigraph returns an igraph object with correct dimensions", {})
