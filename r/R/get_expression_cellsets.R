#' Get expression genes for pathway analysis
#'
#' @param req request parameters
#'   e.g. list(
#'       body = list(
#'         list(geneName = 'Gene1', comparisonType = 'greaterThan', thresholdValue = 0.5),
#'         list(geneName = 'Gene2', comparisonType = 'lessThan', thresholdValue = 0.5),
#'         ...
#'       )
#'   )
#' @param data SeuratObject
#'
#' @return vector of cell ids.
#' @export
#'



getExpressionCellSets <- function(req, data) {
  new_cell_set_data <- getExpressionCellSetIDs(req$body$genesConfig, data)
  keep_ids <- unname(new_cell_set_data$keep_ids)
  cell_set_name <- new_cell_set_data$cell_set_name
  new_cell_set <- list(key = uuid::UUIDgenerate(use.time = TRUE), name = cell_set_name, rootNode = FALSE, color = sample(data@misc$color_pool, 1), cellIds = keep_ids)
  cell_set_class_key <- "scratchpad"
  config <- req$body$config

  insert_set_child_through_api(new_cell_set, config$apiUrl, config$experimentId, cell_set_class_key, config$authJwt)
  return(new_cell_set)
}

getExpressionCellSetIDs <- function(filters, data) {
  # get entrez id for each requested gene name
  gene_names <- sapply(filters, `[[`, "geneName")
  gene_annotations <- data@misc$gene_annotations
  name.match <- match(gene_names, gene_annotations$name)

  # fail if any requested gene names are missing (can't return requested cellset)
  if (anyNA(name.match)) {
    stop("Requested ExpressionCellset with gene name(s) that are not present.")
  }

  enids <- gene_annotations$input[name.match]

  # get expression matrix
  expression_mat <- data[["RNA"]]@data

  # subset cells for each filter
  keep.cells <- rep(TRUE, ncol(data))
  comparisons <- list(greaterThan = `>`, lessThan = `<`)
  comparison_strings <- list(greaterThan = ">", lessThan = "<")

  cell_set_name_vector <- list()
  for (i in seq_along(filters)) {
    filter <- filters[[i]]
    enid <- enids[i]

    # using comparison as functions e.g. `<`(x, y)
    comparison <- comparisons[[filter$comparisonType]]
    pass.filter <- comparison(expression_mat[enid, ], filter$thresholdValue)

    cell_set_name_vector <- c(
      cell_set_name_vector,
      paste0(
        filters[[i]]$geneName,
        comparison_strings[filter$comparisonType],
        filter$thresholdValue
      )
    )
    # keep cells that pass previous filter(s) and current
    keep.cells <- keep.cells & pass.filter
  }

  cell_set_name <- paste(cell_set_name_vector, collapse = ", ")
  new_cell_set_data <- list(keep_ids = data$cells_id[keep.cells], cell_set_name = cell_set_name)
  return(new_cell_set_data)
}

