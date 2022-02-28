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
    pbulk[["RNA"]]@counts,
    samples = data.frame(group),
    genes = pbulk@misc$gene_annotations[, "name", drop = FALSE]
  )


  keep <- edgeR::filterByExpr(y, group = group)

  # get filtered genes to add to results
  disc <- row.names(y)[!keep]

  y <- y[keep, ]
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- gsub("^group", "", colnames(design))

  # silences warning when no replication (1 vs 1 sample comparisons)
  v <- suppressWarnings(limma::voomWithQualityWeights(y, design))
  fit <- limma::lmFit(v)

  # make sure we get the right contrast
  contrast <- "base-background"
  contrast_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)

  # calculate logFC or run DE based on sample size
  nsample <- ncol(y)
  if (nsample == 2) {
    # only logFC if two samples
    res <- data.frame(
      logFC = fit$coefficients[, contrast],
      AveExpr = fit$Amean
    )

    res <- res[order(abs(res$logFC), decreasing = TRUE), ]
  } else {
    # differential expression if 3+ samples
    eb_fit <- limma::eBayes(fit, robust = TRUE)
    res <- limma::topTable(eb_fit, coef = contrast, sort.by = "p", n = Inf)

    # rename columns to match up with wilcoxauc
    res <- res[, c('P.Value', 'logFC', 'AveExpr', 'adj.P.Val')]
    colnames(res) <- c('p_val', 'logFC', 'AveExpr', 'p_val_adj')
  }

  res[disc, ] <- NA
  return(res)
}
