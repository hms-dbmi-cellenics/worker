#' Get background genes for pathway analysis
#'
#' @param req request parameters
#' @param data SeuratObject
#'
#' @return list with gene symbols that are expressed in cells being compared.
#' @export
#'
getBackgroundExpressedGenes <- function(req, data) {

    # minimum total count for gene to be considered expressed
    # from edgeR::filterByExpr
    min.total.count = 15

    # add comparison group to 'custom' slot
    data <- addComparisonGroup(req, data)

    # subset to compared cells
    data <- data[, !is.na(data$custom)]

    # get genes with > min.total.count
    ntot <- Matrix::rowSums(data[['RNA']]@counts)
    keep <- ntot > min.total.count

    genes <- data@misc$gene_annotations$original_name[keep]

    return(list(genes = genes))
}

# adds 'custom' slot to identitify background and base cells
addComparisonGroup <- function(req, data) {
    cells_id <- data$cells_id

    # Remove filtered cells
    background <- intersect(req$body$backgroundCells, cells_id)
    base <- intersect(req$body$baseCells, cells_id)

    data$custom <- NA
    data$custom[cells_id %in% background] <- 'background'
    data$custom[cells_id %in% base] <- 'base'
    return(data)
}
