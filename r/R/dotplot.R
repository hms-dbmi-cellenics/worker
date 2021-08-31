runDotPlot <- function(req, data) {
  type <- req$body$type

  # construct clusters from cell sets
  data$custom <- NA
  cells_id <- data$cells_id
  cell_sets <- req$body$cellSets$children

  for (i in seq_along(cell_sets)) {
    set <- cell_sets[[i]]
    filtered_cells <- intersect(set$cellIds, cells_id)
    data$custom[cells_id %in% filtered_cells] <- i
  }

  if(type=="custom"){
    df <- data@misc$gene_annotations
    genes_subset <- subset(df, toupper(df$name) %in% toupper(req$body$genes))
    genes_subset <- genes_subset[, c("input", "name")]
  }else{
    #type == marker
    #Need the marker features
  }
  dotplot_data <- Seurat::DotPlot(data, features = genes_subset$input, group.by = "custom")$data

  dotplot_data$name <- genes_subset[match(dotplot_data$features.plot, genes_subset$input), "name"]

  dotplot_data <- dotplot_data[, c('avg.exp', 'pct.exp', 'name', 'id')]
  colnames(dotplot_data) <- c('avgExpression', 'cellsPercentage', 'geneName', 'cellSet')

  res <- purrr::transpose(dotplot_data)

  return(res)
}
