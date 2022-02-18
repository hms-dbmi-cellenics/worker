
# Function to retrieve all the doublet score for the current experiment.
# The doublet scores were computing in the data-ingest script. To compute then, we
# use the package Scrublet [1]
#' @export
getDoubletScore <- function(req, data) {
  result <- formatMetadataResult(data, "doublet_scores")
  return(result)
}

# Function to retrieve all the mt-content score for the current experiment.
# The MT-content was computing in the data-ingest script. To compute then, we
# use the function PercentageFeatureSet fom Seurat package. We have been able to identify
# the MT-genes only in Mus musculus, Homo sapiens  or drosophila by grepping "MT[-:]"
#' @export
getMitochondrialContent <- function(req, data) {
  result <- formatMetadataResult(data, "percent.mt")
  return(result)
}

formatMetadataResult <- function(data, column) {

  # check if the experiment has specified column
  if (!column %in% colnames(data@meta.data)) {
      stop(
        generateErrorMessage(
            "R_WORKER_COLUMN_NOT_FOUND",
            paste(column, "is not computed for this experiment.")
        )
      )
  }

  # get the specified column, ordering by cells_id
  cells_id_order <- order(data$cells_id, decreasing = FALSE)
  result <- data@meta.data[cells_id_order, column]
  result <- as.data.frame(result)
  result$cells_id <- data@meta.data$cells_id[cells_id_order]
  result <- result %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    dplyr::select(-cells_id)
  result <- t(unname(as.data.frame(result)))
  result <- purrr::map(result, function(x) {
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


# [1] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell
# Transcriptomic Data. Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3.
# PMID: 30954476; PMCID: PMC6625319.
