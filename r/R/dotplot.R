runDotPlot <- function(req, data) {
  markerGenes <- req$body$markerGenes
  # construct clusters from cell sets
  data$custom <- NA
  cells_id <- data$cells_id
  cell_sets <- req$body$cellSets$children
  subsetCellSets <- req$body$subsetCellSets

  if("cellIds" %in% names(subsetCellSets)){
    subsetIds <- subsetCellSets$cellIds
  }else{
    subsetIds <- list()
    for (i in seq_along(subsetCellSets)) {
      set <- cell_sets[[i]]
      subsetIds <- append(subsetIds,set$cellIds)
    }
  }

  if(length(subsetIds)>0){
    meta_data_subset <- data@meta.data[match(subsetIds, data@meta.data$cells_id), ]
    current_cells <- rownames(meta_data_subset)
    data <- subset(data, cells = current_cells)
  }else{
    return(list())
  }

  for (i in seq_along(cell_sets)) {
    set <- cell_sets[[i]]
    filtered_cells <- intersect(set$cellIds, cells_id)
    data$custom[cells_id %in% filtered_cells] <- set$name
  }

  if (markerGenes) {
      # type == marker
      # Need the marker features
  } else {
      req_genes <- req$body$genes
      annot <- data@misc$gene_annotations
      annot_subset <- subset(annot, toupper(name) %in% toupper(req_genes))
      annot_subset <- annot_subset[, c("input", "name")]
  }

  dotplot_data <- Seurat::DotPlot(data, features = annot_subset$input, group.by = "custom")$data
  #features.plot has the ensemble ids
  dotplot_data$name <- annot_subset[dotplot_data$features.plot, "name"]

  dotplot_data <- dotplot_data %>% transmute(cellSets = as.character(id), geneName=as.character(name), avgExpression = avg.exp, cellsPercentage = pct.exp )

  res <- purrr::transpose(dotplot_data)

  return(res)
}
