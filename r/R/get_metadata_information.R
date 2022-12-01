#' Extract doublet scores
#'
#' @param req list request
#' @param data seurat object
#'
#' @export
getDoubletScore <- function(req, data) {
  result <- formatMetadataResult(data, "doublet_scores")
  return(result)
}


#' Retrieve all the mt-content score for every cell
#'
#' The MT-content was calculated in GEM2S.
#'
#' @param req list
#' @param data Seurat object
#'
#' @export
getMitochondrialContent <- function(req, data) {
  result <- formatMetadataResult(data, "percent.mt")
  return(result)
}


#' Retrieve the number of features for every cell
#'
#' @param req list
#' @param data Seurat object
#'
#' @export
getNGenes <- function(req, data) {
  result <- formatMetadataResult(data, "nFeature_RNA")
  return(result)
}


#' Retrieve the number of UMIs for every cell
#'
#' @param req list
#' @param data Seurat object
#'
#' @export
getNUmis <- function(req, data) {
  result <- formatMetadataResult(data, "nCount_RNA")
  return(result)
}


formatMetadataResult <- function(data, column) {

  # check if the experiment has specified column
  if (!column %in% colnames(data@meta.data)) {
    stop(
      generateErrorMessage(
        error_codes$COLUMN_NOT_FOUND,
        paste(column, "is not computed for this experiment.")
      )
    )
  }

  # create correct size vector with NAs, add values ordered by cell_id
  complete_values <- rep(NA_real_, max(data@meta.data$cells_id) + 1)
  complete_values[data$cells_id + 1] <- data@meta.data[, column]

  # convert to list, replacing NAs with NULLs
  result <- lapply(complete_values, function(x) {
    if (is.na(x)) {
      return(NULL)
    } else {
      return(x)
    }
  })

  # Be aware of possible NA values
  if (any(is.null(result))) {
    warning("There are missing values in the ", column, " results")
  }

  return(result)
}
