#' Title
#'
#' @param req {body: {
#'               markerGenes: True/False determines whether to use marker genes or predefined genes
#'               cellSets: Cellsets to show in the plot. Determines whether to show Louvain/Samples/Custom
#'               subsetCellSets: Cellsets to subset the experiment with.
#'               cellSetsIsAll: Bool value. If true, the experiment will be subsetted to all the cellSets ids.
#'            }
#'            }
#' @param data
#'
#' @return
#' @export
#'
#' @examples
runDotPlot <- function(req, data) {
  markerGenes <- req$body$markerGenes
  data$custom <- NA
  cell_sets <- req$body$cellSets$children
  subsetCellSets <- req$body$subsetCellSets
  cell_sets_is_all <- req$body$cellSetsIsAll

  if (length(cell_sets) < 1) {
    message("The requested Cell Sets are empty. Returning empty results.")
    return(list())
  }

  #Construct ids to subset object
  if (!cell_sets_is_all) {
    subsetIds <- subsetCellSets$cellIds
  } else {
    subsetIds <- list()
    for (i in seq_along(subsetCellSets$children)) {
      set <- cell_sets[[i]]
      subsetIds <- append(subsetIds, set$cellIds)
    }
  }

  #Subset seurat object
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
  for (i in seq_along(cell_sets)) {
    set <- cell_sets[[i]]
    filtered_cells <- intersect(set$cellIds, cells_id)
    if (set$name %in% data$custom) {
      data$custom[cells_id %in% filtered_cells] <- paste0(set$name, i)
    } else {
      data$custom[cells_id %in% filtered_cells] <- set$name
    }
  }
  # If NA values are left in the group, dotplot function will fail.
  data <- subset(data, subset = custom != "NA")

  #Get marker genes or requested gene names.
  if (markerGenes) {
    nFeatures <- req$body$nGenes
    all_markers <- getTopMarkerGenes(nFeatures, data, cell_sets)
    features <- as.data.frame(getMarkerNames(data, all_markers))
    rownames(features) <- features$input
  } else {
    req_genes <- req$body$genes
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
