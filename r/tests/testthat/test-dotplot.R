mock_req <- function() {
    req <- list(
        body = list(
            useMarkerGenes = TRUE,
            numberOfMarkers = 5,
            applyFilter = FALSE,
            filter_by = list(
                cellIds = 0:100
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

test_that("useMarkerGenes works", {
    data <- mock_scdata()
    req <- mock_req()

    runDotPlot(req, data)
})
