#' Cell Cycle Scoring Function
#'
#' This function performs cell cycle scoring on scRNA-seq data. Uses the Seurat method as default.
#'
#' @param req The request object contains API URL, experiment ID, authentication JWT.
#' @param scdata SeuratObject
#'
#' @return Formatted cell sets based on cell cycle scoring.
#'
#' @export
#'
cellCycleScoring <- function(req, scdata) {
  cellSets <- run_cell_cycle_scoring(scdata)

  message("formatting cellsets")

  formatted_cell_sets <-
    format_phase_cellsets(cellSets, scdata@misc$color_pool)
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

#' run_cell_cycle_scoring
#'
#' Uses the Seurat method to calculate cell cycle scores.
#'
#' @param scdata SeuratObject
#'
#' @return data.frame with cell ids and cycle stage prediction, or undetermined if there is an error.
#' @export
#'
run_cell_cycle_scoring <- function(scdata) {
  message("Running Cell Cycle Scoring")

  s_gene_names <- Seurat::cc.genes$s.genes
  g2m_gene_names <- Seurat::cc.genes$g2m.genes

  s_gene_ids <-
    na.omit(scdata@misc$gene_annotations[match(tolower(s_gene_names), tolower(scdata@misc$gene_annotations$name)), "input"])
  g2m_gene_ids <-
    na.omit(scdata@misc$gene_annotations[match(tolower(g2m_gene_names), tolower(scdata@misc$gene_annotations$name)), "input"])

  s_genes <- c(s_gene_names,s_gene_ids)
  g2m_genes <- c(g2m_gene_names,g2m_gene_ids)

  phase <- tryCatch({
    Seurat::CellCycleScoring(
      scdata,
      s.features = s_genes,
      g2m.features = g2m_genes,
      set.ident = TRUE
    )@meta.data$Phase
  }, error = function(err) {
    rep("Undetermined", ncol(scdata))
  })

  cellSets <-
    data.frame(cluster = phase,
               cell_ids = scdata@meta.data$cells_id)

  return(cellSets)
}

#' Format cell cycle phase cellsets
#'
#' @param cell_sets data.frame with cell ids and phase
#' @param color_pool list
#' @param name cell cycle phase name
#'
#' @return list of structured cellsets
#' @export
#'
format_phase_cellsets <- function(cell_sets,
                                  color_pool) {
  message("Formatting cluster cellsets.")

  cell_sets_object <-
    list(
      key = "Phase",
      name = "Phase",
      rootNode = TRUE,
      type = "cellSets",
      children = list()
    )
  for (cluster in unique(cell_sets$cluster)) {
    cells <- cell_sets[cell_sets$cluster == cluster, "cell_ids"]

    new_set <- list(
      key = paste0("Phase", "-", cluster),
      name = cluster,
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
