#' Extract normalized expression matrix
#'
#' This function subsets the Seurat object using the cell IDs received from the
#' python worker, then extracts and returns the normalized expression matrix.
#'
#' Data can be subsetted in the python worker according to one of clusters, samples,
#' and metadata groups, or a combination of them.
#'
#' @param req { body: { subsetBy: cell ids to subset the matrix with } }
#' @param data SeuratObject
#'
#' @return
#' @export
#'
GetNormalizedExpression <- function(req, data) {
  subset_ids <- req$body$subsetBy
  apply_subset <- req$body$applySubset

  if (apply_subset && length(subset_ids) == 0) {
    stop(
      generateErrorMessage(
        error_codes$EMPTY_CELL_SET,
        "No cells match requested filters."
      )
    )
  }

  message("Extracting normalized expression matrix")

  if (apply_subset) {
    message("Number of cells before subsetting: ", ncol(data))
    data <- subsetIds(data, cells_id = subset_ids)
  } else {
    message("No subsetting specified, sending the whole matrix")
  }

  matrix <- as.data.frame(Seurat::GetAssayData(data, slot = "data", assay = "RNA"))

  message("Number of cells in matrix to return: ", ncol(matrix))

  matrix <- tibble::rownames_to_column(matrix, var = " ")

  # vroom_write compresses the file based on filename extension
  vroom::vroom_write(matrix, TMP_RESULTS_PATH_GZ, delim = ",", quote = "none")

  return(TMP_RESULTS_PATH_GZ)
}
