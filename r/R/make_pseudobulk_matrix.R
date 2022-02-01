#' Make aggregated matrix for pseudo bulk differential expression
#'
#' The SeuratObject requires a `custom` meta.data slot, specifying which cells
#' belong to the `base` and `background` groups. It also needs a `samples`
#' meta.data slot.
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

  rownames(agg) <- rownames(counts)
  colnames(agg) <- levels(pbulk_groups)

  # recover metadata
  pbulk_metadata <- colnames(agg)
  pbulk_metadata <- data.table::tstrsplit(pbulk_metadata, split = "_")
  samples <- pbulk_metadata[[1]]
  custom <- pbulk_metadata[[2]]

  # create seurat, and add metadata
  pbulk <- Seurat::CreateSeuratObject(agg)
  pbulk$samples <- samples
  pbulk$custom <- custom
  pbulk@misc$gene_annotations <- gene_annotations

  return(pbulk)
}
