mock_req <- function(){
    cellSets = list(children=list(louvain1=list(cellIds=0:39),louvain2=list(cellIds=40:79)))
    req <- list(body=list(nGenes=5,cellSets=cellSets))
    return(req)
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

test_that("Marker Heatmap returns appropiate format", {
    data <- mock_scdata()
    req <- mock_req()

    res <- runMarkerHeatmap(req, data)

    expect_equal(names(res),c("rawExpression","truncatedExpression"))

    # number of rows is number of cells
    expect_equal(ncol(data),nrow(res$rawExpression))

    # returning only at most limit number of genes
    expect_lte(ncol(res$rawExpression),req$body$nGenes*length(req$body$cellSets$children))
})


req <- mock_req()
data <- mock_scdata()
