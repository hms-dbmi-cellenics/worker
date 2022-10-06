#' Subset a Seurat object with the cell ids
#'
#' @param scdata Seurat object
#' @param cells_id Integer vector of cell IDs to keep
#'
#' @return Subsetted Seurat object
#' @export
#'
subsetIds <- function(scdata, cells_id) {
  meta_data_subset <-
    scdata@meta.data[match(cells_id, scdata@meta.data$cells_id), ]
  current_cells <- rownames(meta_data_subset)
  scdata <- subset(scdata, cells = current_cells)
  return(scdata)
}


#' Send cell set to the API
#'
#' This sends a single, new cell set to the API for patching to the cell sets file.
#'
#' @param new_cell_set named list of cell sets
#' @param api_url string URL of the API
#' @param experiment_id string experiment ID
#' @param cell_set_key string cell set UUID
#' @param auth_JWT string authorization token
#'
#' @export
#'
sendCellsetToApi <-
  function(new_cell_set,
           api_url,
           experiment_id,
           cell_set_key,
           auth_JWT) {
    httr_query <- paste0('$[?(@.key == "', cell_set_key, '")]')
    
    new_cell_set$cellIds <- as.list(new_cell_set$cellIds)

    children <-
      list(list("$insert" = list(index = "-", value = new_cell_set)))

    httr::PATCH(
      paste0(api_url, "/v2/experiments/", experiment_id, "/cellSets"),
      body = list(list(
        "$match" = list(query = httr_query,
                        value = list("children" = children))
      )),
      encode = "json",
      httr::add_headers("Content-Type" = "application/boschni-json-merger+json",
                        "Authorization" = auth_JWT)
    )
  }

#' Update the cell sets through the API
#'
#' Used when re-clustering, cell sets are replaced.
#'
#' @param cell_sets_object
#' @param api_url
#' @param experiment_id
#' @param cell_set_key
#' @param auth_JWT
#'
#' @return
#' @export
#'
updateCellSetsThroughApi <-
  function(cell_sets_object,
           api_url,
           experiment_id,
           cell_set_key,
           auth_JWT) {
    httr_query <- paste0("$[?(@.key == \"", cell_set_key, "\")]")

    httr::PATCH(
      paste0(api_url, "/v2/experiments/", experiment_id, "/cellSets"),
      body = list(list(
        "$match" = list(query = httr_query, value = list("$remove" = TRUE))
      ),
      list("$prepend" = cell_sets_object)),
      encode = "json",
      httr::add_headers("Content-Type" = "application/boschni-json-merger+json",
                        "Authorization" = auth_JWT)
    )
  }
