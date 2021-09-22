runDotPlot <- function(req, data) {
  type <- req$body$type
  req_genes <- req$body$genes

  # construct clusters from cell sets
  data$custom <- NA
  cells_id <- data$cells_id
  cell_sets <- req$body$cellSets$children

  for (i in seq_along(cell_sets)) {
    set <- cell_sets[[i]]
    filtered_cells <- intersect(set$cellIds, cells_id)
    data$custom[cells_id %in% filtered_cells] <- i
  }

  if (type == "custom") {
    annot <- data@misc$gene_annotations
    annot_subset <- subset(annot, toupper(name) %in% toupper(req_genes))
    annot_subset <- annot_subset[, c("input", "name")]
  } else {
    # type == marker
    # Need the marker features
  }

  dotplot_data <- Seurat::DotPlot(data, features = annot_subset$input, group.by = "custom")$data
  #features.plot has the ensemble ids
  dotplot_data$name <- annot_subset[dotplot_data$features.plot, "name"]

  dotplot_data <- dotplot_data[, c("avg.exp", "pct.exp", "name", "id")]
  colnames(dotplot_data) <- c("avgExpression", "cellsPercentage", "geneName", "cellSet")

  res <- purrr::transpose(dotplot_data)

  return(res)
}
