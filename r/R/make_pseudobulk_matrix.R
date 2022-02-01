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
  scdata <- scdata[, !is.na(scdata@meta.data$custom)]

  counts <- scdata[["RNA"]]@counts
  gene_annotations <- scdata@misc$gene_annotations

  groups <- factor(scdata@meta.data$samples)
  group_key <- getSampleGroupKey(scdata)

  agg <- presto::sumGroups(counts, groups, MARGIN = 1)
  agg <- Matrix::Matrix(agg, sparse = TRUE)
  agg <- Matrix::t(agg)

  #row/colnames are lost in aggregation
  rownames(agg) <- rownames(counts)
  colnames(agg) <- levels(groups)

  # recover original metadata
  metadata <- data.frame(samples = colnames(agg))
  # makes sure that original custom value goes to correct sample
  metadata <- dplyr::inner_join(metadata, group_key, by = "samples")
  rownames(metadata) <- metadata$samples

  # create seurat, and add metadata
  pbulk <- Seurat::CreateSeuratObject(agg)
  pbulk <- AddMetaData(pbulk, metadata)
  pbulk@misc$gene_annotations <- gene_annotations

  return(pbulk)
}



#' create sample group key table
#'
#' Makes a table with the original sample and custom values, allowing
#' to not have to trust in the groupings/sample orders.
#'
#' @param scdata A filtered SeuratObject
#'
#' @return a data.frame with a row per sample/custom
#' @export
#'
#' @examples
getSampleGroupKey <- function(scdata) {
  metadata <- unique(dplyr::select(scdata@meta.data, samples, custom))
  rownames(metadata) <- NULL

  return(metadata)
}
