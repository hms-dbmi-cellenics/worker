mock_gene_results <- function() {
    data(exprs, package = 'presto')
    data(y, package = 'presto')

    result <- presto::wilcoxauc(exprs, y)
    result <- result[result$group == "A", ]
    rownames(result) <- result$gene_names <- result$feature
    result <- result[, c("gene_names", "pval", "logFC", "pct_in", "pct_out", "padj", "auc")]
    colnames(result) <- list("gene_names", "p_val", "logFC", "pct_1", "pct_2", "p_val_adj", "auc")
    return(result)
}


test_that("applyFilters works with gene with exact match", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = 'gene_names', expression = 'G5'))

  gene_results <- applyFilters(gene_results, filters)

  expect_equal(gene_results$gene_names, 'G5')
})

test_that("applyFilters works with gene with case insensitive match", {
    gene_results <- mock_gene_results()
    filters <- list(list(columnName = 'gene_names', expression = 'g5'))

    gene_results <- applyFilters(gene_results, filters)

    expect_equal(gene_results$gene_names, 'G5')
})

test_that("applyFilters returns all results that match gene search pattern", {
    gene_results <- mock_gene_results()
    filters <- list(list(columnName = 'gene_names', expression = 'G'))

    filt_results <- applyFilters(gene_results, filters)

    expect_equal(filt_results, gene_results)
})


test_that("applyFilters works with single numeric filter", {
    gene_results <- mock_gene_results()
    filter_up <- list(list(columnName = 'logFC', comparison = 'greaterThan', value = 0))
    filter_dn <- list(list(columnName = 'logFC', comparison = 'lessThan', value = 0))

    results_up <- applyFilters(gene_results, filter_up)
    results_dn <- applyFilters(gene_results, filter_dn)

    expect_equal(results_up, gene_results[gene_results$logFC > 0, ])
    expect_equal(results_dn, gene_results[gene_results$logFC < 0, ])
})

test_that("running applyFilters sequentially is equivalent to running together", {
    gene_results <- mock_gene_results()
    filter_gene <- list(columnName = 'gene_names', expression = 'G1')
    filter_up <- list(columnName = 'logFC', comparison = 'greaterThan', value = 0)

    # apply sequentially
    results_seq <- applyFilters(gene_results, list(filter_gene))
    results_seq <- applyFilters(results_seq, list(filter_up))

    # apply together
    results_together <- applyFilters(gene_results, list(filter_gene, filter_up))

    expect_equal(results_seq, results_together)
})

