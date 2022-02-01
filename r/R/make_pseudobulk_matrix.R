#' Make aggregated matrix for pseudo bulk differential expression
#'
#'
#' @param scdata a SeuratObject,
#'
#' @return a SeuratObject with counts aggregated by sample for one cluster
#' @export
#'
makePseudobulkMatrix <- function(scdata) {
  counts <- scdata[["RNA"]]@counts
  gene_annotations <- scdata@misc$gene_annotations

  # create groups for aggregation
  pbulk_groups <-
    factor(paste(scdata$samples, scdata$custom, sep = "_"))

  agg <- presto::sumGroups(counts, pbulk_groups, MARGIN = 1)
  agg <- Matrix::Matrix(agg, sparse = TRUE)
  agg <- Matrix::t(agg)

  row.names(agg) <- row.names(counts)
  colnames(agg) <- levels(pbulk_groups)

  # recover metadata
  pbulk_groups <- colnames(agg)
  pbulk_groups <- transpose(strsplit(pbulk_groups, split = "_"))
  samples <- pbulk_groups[[1]]
  custom <- pbulk_groups[[2]]

  # create seurat, and add metadata
  pbulk <- CreateSeuratObject(agg)
  pbulk$samples <- samples
  pbulk$custom <- custom
  pbulk@misc$gene_annotations <- gene_annotations

  pbulk

}
