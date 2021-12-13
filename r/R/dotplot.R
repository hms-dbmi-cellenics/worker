#' Run Request for Dot Plot
#'
#' @param req {body: {
#'               useMarkerGenes: TRUE/FALSE determines whether to use marker genes or predefined genes
#'               numberOfMarkers: Int. Number of marker genes to use
#'               customGenesList: List of Strings. List of marker genes to use
#'               groupBy: Cellsets to show in the plot. Determines whether to show Louvain/Samples/Custom
#'               filterBy: Cellsets to subset the experiment with.
#'               isFilterByAll: Bool value. If true, the experiment will be subsetted to all the cellSets ids.
#'              }
#'            }
#' @param data SeuratObject
#'
#' @return
#' @export
#'
#' @examples
runDotPlot <- function(req, data) {
  use_marker_genes <- req$body$useMarkerGenes
  group_by_cell_sets <- req$body$groupBy$children
  filter_by <- req$body$filterBy
  apply_filter <- req$body$applyFilter

  if (length(group_by_cell_sets) < 1) {
    message("The requested Cell Sets are empty. Returning empty results.")
    return(list())
  }

  # subset object to requested cells
  if (apply_filter) {
    subset_ids <- filter_by$cellIds

    if (!length(subset_ids)) {
      message("The ids to subset the object are empty. Returning empty results.")
      return(list())
    }

    data <- subset_ids(data, subset_ids)
  }

  # Construct the dotplot_groups slot
  data$dotplot_groups <- NA

  # This covers a border case where two cell_sets have the same name (but different ID). Can happen in scratchpad
  cell_set_names <- make.unique(sapply(group_by_cell_sets, `[[`, "name"))

  for (i in seq_along(group_by_cell_sets)) {
    cell_set <- group_by_cell_sets[[i]]
    cell_set_name <- cell_set_names[i]
    filtered_cells <- intersect(cell_set$cellIds, cells_id)
    data$dotplot_groups[cells_id %in% filtered_cells] <- cell_set_name
  }

  # If NA values are left in the group, dotplot function will fail.
  subset_cells <- colnames(data)[!is.na(data$dotplot_groups)]
  data <- subset(data, cells = subset_cells)

  # Get marker genes or requested gene names.
  if (use_marker_genes) {
    num_features <- req$body$numberOfMarkers
    all_markers <- getTopMarkerGenes(num_features, data, group_by_cell_sets)
    features <- as.data.frame(getMarkerNames(data, all_markers))
    rownames(features) <- features$input

  } else {
    req_genes <- req$body$customGenesList
    annot <- data@misc$gene_annotations
    annot_subset <- subset(annot, toupper(name) %in% toupper(req_genes))
    features <- annot_subset[, c("input", "name")]
  }

  dotplot_data <- Seurat::DotPlot(data, features = features$input, group.by = "dotplot_groups")$data
  # features.plot has the ensemble ids: get gene symbols
  dotplot_data$name <- features[dotplot_data$features.plot, "name"]
  dotplot_data <- dotplot_data[stringr::str_order(dotplot_data$id, numeric = TRUE), ]
  dotplot_data <- dotplot_data %>%
    transmute(cellSets = as.character(id),
              geneName = as.character(name),
              avgExpression = avg.exp.scaled,
              cellsPercentage = pct.exp)

  res <- purrr::transpose(dotplot_data)
  return(res)
}
