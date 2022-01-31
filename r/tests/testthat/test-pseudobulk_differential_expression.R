mock_pbulk <- function() {

    # mock pseudobulk SeuratObject with annotation
    pbmc_raw <- read.table(
        file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
        as.is = TRUE
    )

    pbmc_raw <- as(pbmc_raw, 'dgCMatrix')

    groups <- rep(letters[1:4], each = 20)
    pbulk <- Matrix.utils::aggregate.Matrix(Matrix::t(pbmc_raw), groupings = groups, fun = "sum")
    pbulk <- Matrix::t(pbulk)

    pbulk <- CreateSeuratObject(counts = pbulk)
    pbulk@misc$gene_annotations <- data.frame(
        input = paste0("ENSG", seq_len(nrow(pbulk))),
        name = row.names(pbulk),
        row.names = paste0("ENSG", seq_len(nrow(pbulk)))
    )

    # add comparison group
    pbulk$custom <- c('base', 'base', 'background', 'background')


    return(pbulk)
}

test_that("runPseudobulkDE runs", {
    pbulk <- mock_pbulk()

    res <- runPseudobulkDE(pbulk)
    expect_equal(class(res), 'data.frame')
})
