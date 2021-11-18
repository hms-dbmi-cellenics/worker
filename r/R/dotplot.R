#' Title
#'
#' @param req {body: {
#'               useMarkerGenes: True/False determines whether to use marker genes or predefined genes
#'               numberOfMarkers: Int. Number of marker genes to use
#'               customGenesList: List of Strings. List of marker genes to use
#'               groupBy: Cellsets to show in the plot. Determines whether to show Louvain/Samples/Custom
#'               filterBy: Cellsets to subset the experiment with.
#'               isFilterByAll: Bool value. If true, the experiment will be subsetted to all the cellSets ids.
#'              }
#'            }
#' @param data
#'
#' @return
#' @export
#'
#' @examples
runDotPlot <- function(req, data) {
  useMarkerGenes <- req$body$useMarkerGenes
  data$custom <- NA
  group_by_cell_sets <- req$body$groupBy$children
  filter_by <- req$body$filterBy
  is_filter_by_all <- req$body$isFilterByAll

  if (length(group_by_cell_sets) < 1) {
    message("The requested Cell Sets are empty. Returning empty results.")
    return(list())
  }

  #Collect ids to subset object
  if (!is_filter_by_all) {
    subsetIds <- filter_by$cellIds
  } else {
    subsetIds <- list()
    for (i in seq_along(filter_by)) {
      filter_parent <- filter_by[[i]]$children
      for(j in seq_along(filter_parent)) {
        subsetIds <- append(subsetIds, filter_parent[[j]]$cellIds)
      }
    }
  }

  #Subset cell Ids
  if (length(subsetIds)) {
    meta_data_subset <- data@meta.data[match(subsetIds, data@meta.data$cells_id), ]
    current_cells <- rownames(meta_data_subset)
    data <- subset(data, cells = current_cells)
    cells_id <- data$cells_id
  } else {
    message("The ids to subset the object are empty. Returning empty results.")
    return(list())
  }

  #Construct the custom slot
  for (i in seq_along(group_by_cell_sets)) {
    cell_set <- group_by_cell_sets[[i]]
    filtered_cells <- intersect(cell_set$cellIds, cells_id)
    if (cell_set$name %in% data$custom) {
      data$custom[cells_id %in% filtered_cells] <- paste0(cell_set$name, i)
    } else {
      data$custom[cells_id %in% filtered_cells] <- cell_set$name
    }
  }
  # If NA values are left in the group, dotplot function will fail.
  data <- subset(data, subset = custom != "NA")

  #Get marker genes or requested gene names.
  if (useMarkerGenes) {
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

  dotplot_data <- Seurat::DotPlot(data, features = features$input, group.by = "custom")$data
  # features.plot has the ensemble ids
  dotplot_data$name <- features[dotplot_data$features.plot, "name"]
  dotplot_data <- dotplot_data[stringr::str_order(dotplot_data$id, numeric = TRUE), ]
  dotplot_data <- dotplot_data %>% transmute(cellSets = as.character(id), geneName = as.character(name), avgExpression = avg.exp, cellsPercentage = pct.exp)

  res <- purrr::transpose(dotplot_data)
  return(res)
}
