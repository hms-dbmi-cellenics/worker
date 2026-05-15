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

  logcounts <- Seurat::GetAssayData(data, layer = "data", assay = "RNA")
  logcounts <- as(logcounts, "dgCMatrix")
  logcounts <- as.data.frame(logcounts)

  message("Number of cells in matrix to return: ", ncol(logcounts))

  logcounts <- tibble::rownames_to_column(logcounts, var = " ")

  # fwrite compresses the file based on filename extension
  data.table::fwrite(
    logcounts,
    TMP_RESULTS_PATH_GZ,
    quote = FALSE,
    row.names = TRUE
  )

  return(TMP_RESULTS_PATH_GZ)
}
