runPseudobulkDE <- function(scdata) {
    library(edgeR)
    group <- scdata$custom

    y <- edgeR::DGEList(
        scdata[['RNA']]@counts,
        samples = data.frame(group),
        genes = scdata@misc$gene_annotations[, "name", drop = FALSE])


    y <- y[edgeR::filterByExpr(y, group=group),]
    y <- edgeR::calcNormFactors(y)

    design <- stats::model.matrix(~0 + group)

    v <- limma::voomWithQualityWeights(y, design)
    fit <- limma::lmFit(v)

    # make sure we get the right contrast
    contrast <- 'base-background'
    contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)

    eb_fit <- limma::eBayes(fit, robust=TRUE)

    res <- limma::topTable(eb_fit, coef = contrast, sort.by="p", n=Inf)
    return(res)
}
