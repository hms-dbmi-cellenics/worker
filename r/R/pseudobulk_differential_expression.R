#' Run Pseudobulk Differential Expression
#'
#' @param pbulk a SeuratObject with counts aggregated by sample for one cluster.
#'
#' @return a data.frame of differential expression results.
#' @export
#'
runPseudobulkDE <- function(pbulk) {
    group <- pbulk$custom

    y <- edgeR::DGEList(
        pbulk[['RNA']]@counts,
        samples = data.frame(group),
        genes = pbulk@misc$gene_annotations[, "name", drop = FALSE])

    keep <- edgeR::filterByExpr(y, group=group)

    # get filtered genes to add to results
    disc <- row.names(y)[!keep]

    y <- y[keep,]
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~0 + group)
    colnames(design) <- gsub('^group', '', colnames(design))

    v <- limma::voomWithQualityWeights(y, design)
    fit <- limma::lmFit(v)

    # make sure we get the right contrast
    contrast <- 'base-background'
    contrast_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)

    eb_fit <- limma::eBayes(fit, robust=TRUE)
    res <- limma::topTable(eb_fit, coef = contrast, sort.by="p", n=Inf)

    # add filtered genes so that searchable
    res[disc, ] <- NA
    res[disc, 'name'] <- pbulk@misc$gene_annotations[disc, 'name']


    return(res)
}
