mock_req <- function() {
    req <- list(
        body = list(
            useMarkerGenes = TRUE,
            customGenesList = c('PF4', 'CST3'),
            numberOfMarkers = 5,
            applyFilter = FALSE,
            filterBy = list(
                cellIds = c(0:25, 51:75)
            ),
            groupBy = list(
                children = list(
                    list(
                        name = 'hello',
                        cellIds = 0:50
                    ),
                    list(
                        name = 'hello2',
                        cellIds = 51:100
                    )
                )
            )
        )
    )
}

mock_scdata <- function() {
    data("pbmc_small", package = 'SeuratObject', envir = environment())
    pbmc_small$cells_id <- 0:(ncol(pbmc_small)-1)
    pbmc_small@misc$gene_annotations <- data.frame(
        input = row.names(pbmc_small),
        name = row.names(pbmc_small),
        row.names = row.names(pbmc_small)
    )
    return(pbmc_small)
}

test_that("dotplot generates the expected list format", {
    data <- mock_scdata()
    req <- mock_req()

    res <- runDotPlot(req, data)
    item <- res[[1]]
    expect_named(item, c('cellSets', 'geneName', 'avgExpression', 'cellsPercentage'))
    expect_type(item$cellSets, 'character')
    expect_type(item$geneName, 'character')
    expect_type(item$avgExpression, 'double')
    expect_type(item$cellsPercentage, 'double')
})

test_that("useMarkerGenes works", {
    data <- mock_scdata()
    req <- mock_req()

    res <- runDotPlot(req, data)
    expect_snapshot(res)
})

test_that("customGenesList is used if useMarkerGenes is FALSE", {
    data <- mock_scdata()
    req <- mock_req()
    req$body$useMarkerGenes <- FALSE

    res <- runDotPlot(req, data)
    genes_used <- unique(sapply(res, `[[`, 'geneName'))

    expect_true(all(genes_used %in% req$body$customGenesList))
})

test_that("subsetting is applied and changes dotpot output", {
    data <- mock_scdata()
    req <- mock_req()

    res_unfilt <- runDotPlot(req, data)

    req$body$applyFilter <- TRUE
    res_filt <- runDotPlot(req, data)

    expect_false(identical(res_unfilt, res_filt))
})
