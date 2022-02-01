

make_pbulk_seurat <- function(scdata) {
  # create pseudobulk groupings
  scdata$pbulk_groups <-
    paste(scdata$groups, scdata$custom, sep = "_")

  # extract annotations,
  gene_annotations <- scdata@misc$gene_annotations

  pbulk <- AggregateExpression(
    scdata,
    slot = "counts",
    group.by = "pbulk_groups",
    return.seurat = TRUE,
    verbose = FALSE
  )

  # add annotations to aggregated object
  pbulk@misc$gene_annotations <- gene_annotations

  # fix meta.data
  pbulk_groups <- rownames(pbulk@meta.data)
  pbulk_groups <- transpose(strsplit(pbulk_groups, split = "_"))

  samples <- pbulk_groups[[1]]
  custom <- pbulk_groups[[2]]

  pbulk@meta.data$samples <- samples
  pbulk@meta.data$custom <- custom

  pbulk
}
