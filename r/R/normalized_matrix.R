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

  if (length(subset_ids) == 0) {
    stop(
      generateErrorMessage(
        error_codes$EMPTY_CELL_SET,
        "No cells match requested filters."
      )
    )
  }

  message("Extracting normalized expression matrix")
  message("Number of cells before subsetting: ", ncol(data))

  if (apply_filter == TRUE) {
    data <- subsetIds(data, cells_id = subset_ids)
  }
  norm_matrix <- as.data.frame(Seurat::GetAssayData(data, slot = "data", assay = "RNA"))

  message("Number of cells after subsetting: ", ncol(norm_matrix))
  message("Extracting normalized expression matrix from whole data")

  return(norm_matrix)
}