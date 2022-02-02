#' Make aggregated matrix for pseudo bulk differential expression
#'
#' The SeuratObject requires a `custom` meta.data slot, specifying which cells
#' belong to the `base` and `background` groups. It filters out cells with `NA` in `custom`.
#' Also, a `samples` meta.data slot is required, specifying samples.
#'
#' @param scdata a SeuratObject,
#'
#' @return a SeuratObject with counts aggregated by sample for one cluster
#' @export
#'
makePseudobulkMatrix <- function(scdata) {
  # filter out cells not in base/background groups
  scdata <- scdata[, !is.na(scdata$custom)]

  counts <- scdata[["RNA"]]@counts
  gene_annotations <- scdata@misc$gene_annotations

  groups <- factor(scdata$samples)
  # aggregate over samples
  samples <- scdata$samples
  samples <- factor(samples, levels = unique(samples))

  agg <- presto::sumGroups(counts, samples, MARGIN = 1)
  agg <- Matrix::Matrix(agg, sparse = TRUE)
  agg <- Matrix::t(agg)

  #row/colnames are lost in aggregation
  rownames(agg) <- rownames(counts)
  colnames(agg) <- levels(samples)

  # recover original metadata
  metadata <- data.frame(samples = colnames(agg),
                         custom = scdata$custom[!duplicated(samples)],
                         row.names = colnames(agg))

  # create seurat, and add metadata
  pbulk <- Seurat::CreateSeuratObject(agg, meta.data = metadata)
  pbulk@misc$gene_annotations <- gene_annotations

  return(pbulk)
