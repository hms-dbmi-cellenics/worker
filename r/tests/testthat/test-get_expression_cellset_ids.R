mock_scdata <- function() {

    pbmc_raw <- read.table(
        file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
        as.is = TRUE
    )

    gene_annotations <- data.frame(
        input = paste0("ENSG", seq_len(nrow(pbmc_raw))),
        name = row.names(pbmc_raw),
        row.names = paste0("ENSG", seq_len(nrow(pbmc_raw)))
    )

    row.names(pbmc_raw) <- gene_annotations$input


    pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)

    pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
    pbmc_small@misc$gene_annotations <- gene_annotations
    return(pbmc_small)
}

test_that("getExpressionCellsetIDs with a single expression filter returns corect cell ids", {
    data <- mock_scdata()
    req <- list(
        body = list(
            list(geneName = 'MS4A1', comparisonType = 'greaterThan', thresholdValue = 0.5)
        )
    )
    res <- getExpressionCellsetIDs(req, data)

    annot <- data@misc$gene_annotations
    enid <- annot$input[annot$name == req$body[[1]]$geneName]
    expected <- data$cells_id[data[['RNA']]@data[enid, ] > 0.5]
    expect_equal(res, expected)
})

test_that("getExpressionCellsetIDs with two expression filters returns corect cell ids", {
    data <- mock_scdata()
    req <- list(
        body = list(
            list(geneName = 'MS4A1', comparisonType = 'greaterThan', thresholdValue = 0.5),
            list(geneName = 'CD79B', comparisonType = 'greaterThan', thresholdValue = 0.5)
        )
    )
    res <- getExpressionCellsetIDs(req, data)

    # gives same results as manual
    annot <- data@misc$gene_annotations
    enid1 <- annot$input[annot$name == req$body[[1]]$geneName]
    enid2 <- annot$input[annot$name == req$body[[2]]$geneName]
    expected <- data$cells_id[data[['RNA']]@data[enid1, ] > 0.5 & data[['RNA']]@data[enid2, ] > 0.5]

    expect_equal(res, expected)

    # order invariant
    req2 <- list(
        body = list(
            list(geneName = 'CD79B', comparisonType = 'greaterThan', thresholdValue = 0.5),
            list(geneName = 'MS4A1', comparisonType = 'greaterThan', thresholdValue = 0.5)
        )
    )

    res2 <- getExpressionCellsetIDs(req2, data)
    expect_equal(res, res2)
})


test_that("getExpressionCellsetIDs returns no cell ids if filters remove everything", {
    data <- mock_scdata()
    req <- list(
        body = list(
            list(geneName = 'MS4A1', comparisonType = 'greaterThan', thresholdValue = 0.5),
            list(geneName = 'MS4A1', comparisonType = 'lessThan', thresholdValue = 0.5)
        )
    )
    res <- getExpressionCellsetIDs(req, data)
    expect_length(res, 0)
})


test_that("getExpressionCellsetIDs returns all cell ids if filters include everything", {
    data <- mock_scdata()
    req <- list(
        body = list(
            list(geneName = 'MS4A1', comparisonType = 'greaterThan', thresholdValue = -1)
        )
    )
    res <- getExpressionCellsetIDs(req, data)
    expect_equal(res, data$cells_id)
})



test_that("getExpressionCellsetIDs fails if requested geneNames are not present in SeuratObject", {
    data <- mock_scdata()
    req <- list(
        body = list(
            list(geneName = 'BLAH', comparisonType = 'greaterThan', thresholdValue = 0)
        )
    )

    expect_error(getExpressionCellsetIDs(req, data), "gene name\\(s\\) that are not present[.]")
})
