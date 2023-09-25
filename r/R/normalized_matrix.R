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

  # Write the matrix to a temporary file
  tmp_file <- tempfile()
  data.table::fwrite(data.table::as.data.table(matrix), tmp_file, sep = ",")

  # Read the contents of the file into a raw vector
  raw_data <- readBin(tmp_file, "raw", file.info(tmp_file)$size)

  # Initialize an empty list to store the encoded chunks
  encoded_data <- list()

  # Define the chunk size
  chunk_size <- 2^20  # 1 MB

  # Open a connection to the temporary file
  con <- file(tmp_file, "rb")

  repeat {
    # Read a chunk from the file into a raw vector
    raw_chunk <- readBin(con, "raw", chunk_size)

    # If the chunk is empty, we've reached the end of the file
    if (length(raw_chunk) == 0) break

    # Encode the chunk and append it to the list
    encoded_data[[length(encoded_data) + 1]] <- base64enc::base64encode(raw_chunk)
  }

  # Close the connection and delete the temporary file
  close(con)
  unlink(tmp_file)

  return(encoded_data)
}
