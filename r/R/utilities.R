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
#' This sends a single, new cell set to the API for patching to the cell sets
#' file.
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
        "$match" = list(
          query = httr_query,
          value = list("children" = children)
        )
      )),
      encode = "json",
      httr::add_headers(
        "Content-Type" = "application/boschni-json-merger+json",
        "Authorization" = auth_JWT
      )
    )
  }

#' Update the cell sets through the API
#'
#' Used when re-clustering, cell sets are replaced.
#'
#' @param cell_sets_object list of cellsets to patch
#' @param api_url character - api endpoint url
#' @param experiment_id character
#' @param cell_set_key character
#' @param auth_JWT character
#'
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
      body = list(
        list(
          "$match" = list(query = httr_query, value = list("$remove" = TRUE))
        ),
        list("$prepend" = cell_sets_object)
      ),
      encode = "json",
      httr::add_headers(
        "Content-Type" = "application/boschni-json-merger+json",
        "Authorization" = auth_JWT
      )
    )
  }


#' Ensure is list in json
#'
#' When sending responses as json, Vectors of length 0 or 1 are converted to
#' null and scalar (respectively) Using as.list fixes this, however, long R
#' lists take a VERY long time to be converted to JSON.
#' This function deals with the problematic cases, leaving vector as a vector
#' when it isnt a problem.
#'
#' @param vector
#'
#' @export
#'
ensure_is_list_in_json <- function(vector) {
  if (length(vector) <= 1) {
    return(as.list(vector))
  } else {
    return(vector)
  }
}


#' Add NAs to fill variables for filtered cell ids
#'
#' This function creates a vector of size max(cell_ids) + 1, with NAs in each
#' index that corresponds to a filtered cell and the corresponding value in the
#' ones that were not. It returns the values ordered by cell id by design.
#'
#' @param variable vector of values to complete
#' @param cell_ids integer vector of filtered cell ids
#'
#' @return NA filled vector, cell_id-complete
#' @export
#'
complete_variable <- function(variable, cell_ids) {
  # create correct size vector with NAs, add values ordered by cell_id
  complete_values <- rep(NA_real_, max(cell_ids) + 1)
  complete_values[cell_ids + 1] <- variable
  return(complete_values)
}


add_clusters <- function(scdata, parsed_cellsets) {
  seurat_clusters <- parsed_cellsets[cellset_type == "cluster", c("name", "cell_id")]
  data.table::setnames(seurat_clusters, c("seurat_clusters", "cells_id"))
  scdata@meta.data <- dplyr::left_join(scdata@meta.data, seurat_clusters, by = "cells_id")

  if("scratchpad" %in% parsed_cellsets[["cellset_type"]]) {
    custom_clusters <- parsed_cellsets[cellset_type == "scratchpad", c("name", "cell_id")]
    # create one column for each scratchpad cluster because one cell can be assigned to more than one scratchpad cluster
    custom_clusters_list <- split(custom_clusters, custom_clusters[["name"]])
    for (i in 1:length(custom_clusters_list)) {
      scratchpad_colname <- paste0("scratchpad-", names(custom_clusters_list)[i])
      data.table::setnames(custom_clusters_list[[i]], c(scratchpad_colname, "cells_id"))
      scdata@meta.data <- dplyr::left_join(scdata@meta.data, custom_clusters_list[[i]], by = "cells_id")
    }
  }

  return(scdata)
}

parse_cellsets <- function(cellsets) {
  # filter out elements with length = 0 (e.g. if scratchpad doesn't exist)
  cellsets <- cellsets[sapply(cellsets, length) > 0]

  dt <- purrr::map2_df(cellsets, names(cellsets), ~ cbind(cellset_type = .y, rrapply::rrapply(.x, how = "bind")))
  dt <- data.table::setDT(dt)
  dt <- dt[, setNames(.(unlist(cellIds)), "cell_id"), by = .(key, name, cellset_type)]

  # change cellset type to more generic names
  dt[cellset_type %in% c("louvain", "leiden"), cellset_type := "cluster"]
  dt[!cellset_type %in% c("cluster", "scratchpad", "sample"), cellset_type := "metadata"]

  return(dt)
}


#' Determine the type of features in the annot data frame
#'
#' Classifies the features file columns into either ensemblIds or symbols
#'
#' @param annot data.frame read from features file
#' @return character vector indicating feature types
#'
#' @export
get_feature_types <- function(annot) {
  is_ens_col1 <- startsWith(annot[[1]], "ENS")
  pct_ens_col1 <- sum(is_ens_col1) / nrow(annot)

  is_ens_col2 <- startsWith(annot[[2]], "ENS")
  pct_ens_col2 <- sum(is_ens_col2) / nrow(annot)

  is_ens <- c(pct_ens_col1, pct_ens_col2) >= 0.5

  # reverse case, sym in first and id in second column
  if (!is_ens[1] && is_ens[2]) {
    return(SYM_IDS)
  }

  # regular cases. sum of booleans returns ints. convert to char to string match
  feature_type <- switch(as.character(sum(is_ens)),
                         "0" = SYM_SYM,
                         "1" = IDS_SYM,
                         "2" = IDS_IDS
  )

  return(feature_type)
}


format_sctype_cell_sets <-
  function(data, species, tissue, color_pool) {
    cell_class_key <- paste0("ScType-", tissue, "-", species)

    cell_class <-
      list(
        key = cell_class_key,
        name = cell_class_key,
        rootNode = TRUE,
        type = "cellSets",
        children = list()
      )

    for (i in 1:length(unique(data@meta.data$customclassif))) {
      cell_set_key <- unique(data@meta.data$customclassif)[[i]]

      new_cell_set <- list(
        key = cell_set_key,
        name = cell_set_key,
        rootNode = FALSE,
        type = "cellSets",
        color = sample(data@misc$color_pool, 1),
        cellIds = data@meta.data[data@meta.data$customclassif == cell_set_key, "cells_id"]
      )
      color_pool <- color_pool[-1]
      cell_class$children <- append(cell_class$children, list(new_cell_set))
    }

    return(cell_class)
  }
