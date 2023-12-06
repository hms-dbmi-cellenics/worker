cellCycleScoring <- function(req,data){
  cellSets <- run_cycle_scoring(scdata)

  message("formatting cellsets")

  formatted_cell_sets <-
    format_cluster_cellsets(cellSets, "Phase", scdata@misc$color_pool)
  message("updating through api")

  updateCellSetsThroughApi(
    formatted_cell_sets,
    req$body$apiUrl,
    req$body$experimentId,
    "Phase",
    req$body$authJwt,
    append = FALSE
  )

  return(formatted_cell_sets)
}

run_cell_cycle_scoring <- function(scdata, req){
  message("Running Cell Cycle Scoring")
  s.genes <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes

  phase <- tryCatch({
    Seurat::CellCycleScoring(scdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)@meta.data$Phase
  }, error = function(err) {
    rep("Undetermined",ncol(scdata))
  })

  cellSets <- data.frame(cluster = phase, cell_ids = scored_scdata@meta.data$cells_id)

  return(cellSets)
}

format_cluster_cellsets <- function(cell_sets,
                                    clustering_method,
                                    color_pool,
                                    name = paste0(clustering_method)) {
  message("Formatting cluster cellsets.")

  # careful with capital l on type for the key.
  cell_sets_object <-
    list(
      key = clustering_method,
      name = name,
      rootNode = TRUE,
      type = "cellSets",
      children = list()
    )
  for (cluster in unique(cell_sets$cluster)) {
    cells <- cell_sets[cell_sets$cluster == cluster, "cell_ids"]
    is.num <- !is.na(as.numeric(cluster))
    set_name <- ifelse(is.num, paste("Cluster", cluster), cluster)

    new_set <- list(
      key = paste0(clustering_method, "-", cluster),
      name = set_name,
      rootNode = FALSE,
      type = "cellSets",
      color = color_pool[1],
      cellIds = ensure_is_list_in_json(unname(cells))
    )
    color_pool <- color_pool[-1]
    cell_sets_object$children <-
      append(cell_sets_object$children, list(new_set))
  }
  return(cell_sets_object)
}

ensure_is_list_in_json <- function(vector) {
  if (length(vector) <= 1) {
    return(as.list(vector))
  } else {
    return(vector)
  }
}

