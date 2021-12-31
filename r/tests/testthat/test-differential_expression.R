mock_req <- function() {
    req <- list(
        body = list(
            backgroundCells = 0:50,
            baseCells = 51:100,
            pagination = list(
                orderBy = 'logFC',
                orderDirection = 'DESC',
                offset = 0,
                limit = 50
            )
        )
    )
}

mock_scdata <- function() {
    data("pbmc_small", package = 'SeuratObject', envir = environment())
    pbmc_small$cells_id <- 0:(ncol(pbmc_small)-1)
    pbmc_small@misc$gene_annotations <- data.frame(
        input = paste0('ENSG', seq_len(nrow(pbmc_small))),
        name = row.names(pbmc_small),
        row.names = paste0('ENSG', seq_len(nrow(pbmc_small)))
    )
    return(pbmc_small)
}

test_that("runDE generates the expected return format", {
    data <- mock_scdata()
    req <- mock_req()

    res <- runDE(req, data)

    # number of genes is number of possible DE rows
    expect_equal(res$full_count, nrow(data))

    # returning only at most limit number of genes
    expect_equal(nrow(res$gene_results), req$body$pagination$limit)

    # ordering is correct
    expect_equal(res$gene_results$logFC, sort(res$gene_results$logFC, decreasing = TRUE))

    # have the correct column names
    expect_columns <- c('p_val', 'logFC', 'pct_1', 'pct_2', 'p_val_adj', 'auc', 'gene_names', 'Gene')
    expect_equal(colnames(res$gene_results), expect_columns)
})

test_that("runDE limit won't return more than available genes", {
    data <- mock_scdata()
    req <- mock_req()

    req$body$pagination$limit <- nrow(data) + 50

    res <- runDE(req, data)
    expect_equal(nrow(res$gene_results), nrow(data))
})

test_that("runDE was able to convert from ensembl ids to gene symbols", {
    data <- mock_scdata()
    req <- mock_req()

    res <- runDE(req, data)
    expect_equal(sum(is.na(res$gene_results$gene_names)), 0)
    expect_true(all(res$gene_results$gene_names %in% data@misc$gene_annotations$name))
    expect_equal(length(unique(res$gene_results$gene_names)), nrow(res$gene_results))

})

test_that("runDE works with gene name filter", {
    data <- mock_scdata()
    req <- mock_req()
    req$body$pagination$filters <-
        list(list(columnName = 'gene_names', expression = 'CST3'))


    res <- runDE(req, data)
    expect_equal(row.names(res$gene_results), 'CST3')
})

test_that("runDE works with numeric filters", {
    data <- mock_scdata()
    req <- mock_req()
    req$body$pagination$filters <-
        list(list(columnName = 'logFC', comparison = 'greaterThan', value = 0))


    res <- runDE(req, data)
    expect_true(all(res$gene_results$logFC > 0))
})

test_that("runDE works when no genes match filters", {
    data <- mock_scdata()
    req <- mock_req()
    req$body$pagination$filters <-
        list(list(columnName = "gene_names", expression = "mythical_gene"))

    res <- runDE(req, data)
    expect_equal(nrow(res$gene_results), 0)

    req$body$pagination$filters <-
        list(list(
            columnName = "logFC",
            comparison = "greaterThan",
            value = Inf
        ))

    res <- runDE(req, data)
    expect_equal(nrow(res$gene_results), 0)
})
