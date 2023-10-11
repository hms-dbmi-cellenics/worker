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

  INTERNAL_RESULTS_PATH <- "/data/rResult"

  print("HOLA3")

  tryCatch({
    write.csv(matrix, INTERNAL_RESULTS_PATH)
  }, warning = function(warning_condition) {
    print("warning_conditionDebug")
    print(warning_condition)
  }, error = function(error_condition) {
    print("error_conditionDebug")
    print(error_condition)
  }, finally={})


  print("HOLA5")

  return(INTERNAL_RESULTS_PATH)
}