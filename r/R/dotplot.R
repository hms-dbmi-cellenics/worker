runDotPlot <- function(req, data) {
  type <- req$body$type
  cellSets <- req$body$cellSets$children
  data$custom <- NA

  object_ids <- data$cells_id
  for (i in seq_along(cellSets)) {
    set <- cellSets[[i]]
    filtered_cells <- intersect(set$cellIds, object_ids)
    data$custom[object_ids %in% filtered_cells] <- i
  }

  if(type=="custom"){
    df <- data@misc$gene_annotations
    genesSubset <- subset(df, toupper(df$name) %in% toupper(req$body$genes))
    genesSubset <- genesSubset[, c("input", "name")]
  }else{
    #type == marker
    #Need the marker features
  }
  dotplot_data <- Seurat::DotPlot(data, features = genesSubset$input, group.by = "custom")$data

  dotplot_data$name <- genesSubset[match(dotplot_data$features.plot, genesSubset$input), "name"]
  rownames(dotplot_data) <- c()
  dotplot_data<-unname(dotplot_data)

  format_results <- function(x) {
    list(avgExpression = as.numeric(x[1]), cellsPercentage = as.numeric(x[2]), geneName = x[5], cellSet = x[4])
  }

  res <- apply(dotplot_data, c(1), format_results)

  return(res)
}
