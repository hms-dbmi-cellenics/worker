mock_scratchpad_cellset_object <- function(n) {
  # ensure cellIds is an int vector. same as when created by getExpressionCellSet
  list(
    key = "cell_set_key",
    name = "cell_set_name",
    rootNode = FALSE,
    color = "color",
    cellIds = as.integer(1:n)
  )
}

mock_cellset_from_python <- function(data) {
  cell_sets <- list(
    "louvain" = list(
      list(
        "key" = "louvain-0", "name" = "Cluster 0", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#77aadd", "cellIds" = unname(data$cells_id[1:ceiling(length(data$cells_id) / 3)])
      ),
      list(
        "key" = "louvain-1", "name" = "Cluster 1", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#77aadd", "cellIds" = unname(data$cells_id[(ceiling(length(data$cells_id) / 3 + 1)):(ceiling(length(data$cells_id) / 3 * 2))])
      ),
      list(
        "key" = "louvain-2", "name" = "Cluster 2", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#ee8866", "cellIds" = unname(data$cells_id[(ceiling(length(data$cells_id) / 3 * 2 + 1)):(length(data$cells_id))])
      )
    ),
    "scratchpad" = list(
      list(
        "key" = "scratchpad-0", "name" = "Custom 0", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#77aadd", "cellIds" = unname(sample(data$cells_id, 5))
      ),
      list(
        "key" = "scratchpad-1", "name" = "Custom 1", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#ee8866", "cellIds" = unname(sample(data$cells_id, 10))
      )
    ),
    "sample" = list(
      list(
        "key" = "a636ec18-4ba3-475b-989d-0a5b2", "name" = "P13 Acute MISC",
        "color" = "#77aadd", "cellIds" = unname(data$cells_id[1:length(data$cells_id) / 2])
      ),
      list(
        "key" = "eec701f2-5762-4b4f-953d-6aba8", "name" = "P13 Convalescent MISC",
        "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
      )
    ),
    "MISC_status" = list(
      list(
        "key" = "MISC_status-Acute", "name" = "Acute", "color" = "#77aadd",
        "cellIds" = unname(data$cells_id[1:length(data$cells_id) / 2])
      ),
      list(
        "key" = "MISC_status-Convalescent", "name" = "Convalescent",
        "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
      )
    )
  )

  return(cell_sets)
}

mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  enids <- paste0("ENSG", seq_len(nrow(pbmc_raw)))
  gene_annotations <- data.frame(
    input = enids,
    name = row.names(pbmc_raw),
    original_name = row.names(pbmc_raw),
    row.names = enids
  )

  row.names(pbmc_raw) <- enids
  pbmc_small <- SeuratObject::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations

  pbmc_small <- Seurat::NormalizeData(pbmc_small, normalization.method = "LogNormalize", verbose = FALSE)
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)

  pbmc_small@misc$color_pool <- mock_color_pool(500)
  return(pbmc_small)
}

with_fake_http(test_that("sendCellSetToApi sends patch request", {
  scr_cellset_object <- mock_scratchpad_cellset_object(5)
  expect_PATCH(sendCellsetToApi(scr_cellset_object, "api_url", "experiment_id", "cell_set_key", "auth"))
}))

mock_color_pool <- function(n) {
  paste0("color_", 1:n)
}

with_fake_http(test_that(
  "sendCellSetToApi sends the cellIds as an array, when there are >1 cells",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(10)
    req <- sendCellsetToApi(
      scr_cellset_object,
      "api_url",
      "experiment_id",
      "cell_set_key",
      "auth"
    )

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    # extract the cellIds slot in the request
    cell_ids <-
      req_body[[1]]$`$match`$value$children[[1]]$`$insert`$value$cellIds

    expect_type(cell_ids, "list")
  }
))


with_fake_http(test_that(
  "sendCellSetToApi sends the cellIds as an array, when there is exactly 1 cell",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(1)
    req <- sendCellsetToApi(
      scr_cellset_object,
      "api_url",
      "experiment_id",
      "cell_set_key",
      "auth"
    )

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    # extract the cellIds slot in the request
    cell_ids <-
      req_body[[1]]$`$match`$value$children[[1]]$`$insert`$value$cellIds

    expect_type(cell_ids, "list")
  }
))


test_that("complete_variable fills a vector with NAs where there are filtered cells", {
  mock_variable <- rnorm(50)
  mock_cell_ids <- seq.int(0, 99, 2)

  res <- complete_variable(mock_variable, mock_cell_ids)

  expect_equal(length(res), max(mock_cell_ids) + 1)
})


test_that("complete_variable returns results ordered by increasing cell id", {
  set.seed(10)
  mock_variable <- rnorm(50)
  mock_cell_ids <- seq.int(0, 99, 2)

  # shuffle cell ids
  mock_cell_ids <- sample(mock_cell_ids, size = length(mock_cell_ids), replace = FALSE)

  # sort the variable by cell id
  cell_id_sorted_variable <- mock_variable[order(mock_cell_ids)]

  res <- complete_variable(mock_variable, mock_cell_ids)

  # remove NAs from result vector
  expect_equal(res[!is.na(res)], cell_id_sorted_variable)
})


test_that("parse_cellsets converts the cellsets object in the expected format", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)
  expected_columns <-
    c(
      "key",
      "name",
      "cellset_type",
      "cell_id"
    )

  parsed_cellsets <- parse_cellsets(cell_sets)

  expect_true(data.table::is.data.table(parsed_cellsets))
  expect_equal(colnames(parsed_cellsets), expected_columns)
  expect_equal(length(cell_sets$louvain), length(unique(parsed_cellsets[cellset_type == "cluster", name])))
  expect_equal(sapply(cell_sets$louvain, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "cluster", name]))
  expect_equal(length(cell_sets$scratchpad), length(unique(parsed_cellsets[cellset_type == "scratchpad", name])))
  expect_equal(sapply(cell_sets$scratchpad, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "scratchpad", name]))
  expect_equal(length(cell_sets$sample), length(unique(parsed_cellsets[cellset_type == "sample", name])))
  expect_equal(sapply(cell_sets$sample, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "sample", name]))
  expect_equal(length(cell_sets$MISC_status), length(unique(parsed_cellsets[cellset_type == "metadata", name])))
  expect_equal(sapply(cell_sets$MISC_status, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "metadata", name]))
})


test_that("add_clusters adds cluster information as metadata columns to the seurat object", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)
  parsed_cellsets <- parse_cellsets(cell_sets)

  data <- add_clusters(data, parsed_cellsets)
  expect_true("seurat_clusters" %in% colnames(data@meta.data))
  expect_equal(all(is.na(data@meta.data$seurat_clusters)), FALSE)

  expect_equal(unique(data@meta.data$seurat_clusters), unique(parsed_cellsets[cellset_type == "cluster", name]))
  expect_true(all.equal(parsed_cellsets[cellset_type == "cluster", c("name", "cell_id")],
    data@meta.data[, c("seurat_clusters", "cells_id")],
    check.attributes = FALSE
  ))

  scratchpad_clusters <- unique(parsed_cellsets[cellset_type == "scratchpad", name])
  for (i in 1:length(scratchpad_clusters)) {
    expected_cells <- which(data@meta.data$cells_id %in% parsed_cellsets[name == scratchpad_clusters[i], cell_id])
    expect_true(all(data@meta.data[expected_cells, paste0("scratchpad-", scratchpad_clusters[i])]))
  }
})


test_that("format_sctype_cell_sets correctly format cellset to be sent to the API", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)
  species <- "human"
  tissue <- "Pancreas"
  active_assay <- "RNA"

  expected_cell_class_names <- c("key", "name", "rootNode", "type", "children")

  scale_data <- get_formatted_data(data, active_assay)
  parsed_cellsets <- parse_cellsets(cell_sets)
  data <- add_clusters(data, parsed_cellsets)
  data[[active_assay]]@scale.data <- scale_data

  data <- suppressWarnings(run_sctype(data, active_assay, tissue, species))

  formatted_cell_class <- format_sctype_cell_sets(data, species, tissue)

  expected_classes <- unique(data@meta.data$customclassif)

  expect_equal(formatted_cell_class$name, paste0("ScType-", tissue, "-", species))
  expect_equal(names(formatted_cell_class), expected_cell_class_names)
  expect_equal(sapply(formatted_cell_class$children, function(x) {
    x$name
  }), expected_classes)

  for (i in 1:length(expected_classes)) {
    expect_equal(
      sapply(formatted_cell_class$children, function(x) {x$cellIds})[[i]],
      data@meta.data[which(data@meta.data$customclassif == expected_classes[i]), "cells_id"]
    )
  }
})


with_fake_http(test_that("updateCellSetsThroughApi sends patch request using append = TRUE", {
  expect_PATCH(
    updateCellSetsThroughApi(list(), "api_url", "experiment_id", "cell_set_key", "auth", append = TRUE)
  )
}))


with_fake_http(test_that("updateCellSetsThroughApi sends patch request using append = FALSE", {
  expect_PATCH(
    updateCellSetsThroughApi(list(), "api_url", "experiment_id", "cell_set_key", "auth", append = FALSE)
  )
}))


