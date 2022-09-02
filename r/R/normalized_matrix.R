#' Extract normalized expression matrix
#'
#' This function subsets the Seurat object using the cell IDs received from the
#' python worker, then extracts and returns the normalized expression matrix.
#'
#' Data can be subsetted in the python worker according to one of clusters, samples,
#' and metadata groups, or a combination of them.
#'
#' @param req {body: {
#'               filterBy: Cellsets to subset the experiment with
#'               applyFilter: TRUE/FALSE determines whether to subset the data
#'              }
#'            }
#' @param data SeuratObject
#'
#' @return
#' @export
#'
GetNormalizedExpression <- function(req, data) {
  apply_filter <- req$body$applyFilter
  subset_ids <- req$body$filterBy

  if(length(subset_ids)==0) {
    stop(
      generateErrorMessage(
        error_codes$EMPTY_CELL_SET,
        "No cells match requested filters."
      )
    )
  }

  if (apply_filter == TRUE) {
    message("Number of cells before subsetting: ", ncol(data))
    data <- subsetIds(data, cells_id = subset_ids)
    message("Extracting normalized expression matrix from subsetted data")
    norm_matrix <- as.data.frame(Seurat::GetAssayData(data, slot = "data", assay = "RNA"))
    message("Number of cells after subsetting: ", ncol(norm_matrix))
  } else {
    message("Extracting normalized expression matrix from whole data")
    norm_matrix <- as.data.frame(Seurat::GetAssayData(data, slot = "data", assay = "RNA"))
  }

  return(norm_matrix)
}
