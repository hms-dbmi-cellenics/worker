library("httptest")

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
      "key" = "louvain",
      "name" = "louvain clusters",
      "rootNode" = TRUE,
      "type" = "cellSets",
      "children" = list(
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
      )
    ),
    "scratchpad" = list(
      "key" = "scratchpad",
      "name" = "Custom cell sets",
      "rootNode" = TRUE,
      "type" = "cellSets",
      "children" = list(
        list(
          "key" = "scratchpad-0", "name" = "Custom 0", "rootNode" = FALSE,
          "type" = "cellSets", "color" = "#77aadd", "cellIds" = unname(sample(data$cells_id, 5))
        ),
        list(
          "key" = "scratchpad-1", "name" = "Custom 1", "rootNode" = FALSE,
          "type" = "cellSets", "color" = "#ee8866", "cellIds" = unname(sample(data$cells_id, 10))
        )
      )
    ),
    "sample" = list(
      "key" = "sample",
      "name" = "Samples",
      "rootNode" = TRUE,
      "type" = "metadataCategorical",
      "children" = list(
        list(
          "key" = "a636ec18-4ba3-475b-989d-0a5b2", "name" = "P13 Acute MISC",
          "color" = "#77aadd", "cellIds" = unname(data$cells_id[1:(length(data$cells_id) / 2)])
        ),
        list(
          "key" = "eec701f2-5762-4b4f-953d-6aba8", "name" = "P13 Convalescent MISC",
          "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
        )
      )
    ),
    "MISC_status" = list(
      "key" = "MISC_status",
      "name" = "MISC_status",
      "rootNode" = TRUE,
      "type" = "metadataCategorical",
      "children" = list(
        list(
          "key" = "MISC_status-Acute", "name" = "Acute", "color" = "#77aadd",
          "cellIds" = unname(data$cells_id[1:(length(data$cells_id) / 2)])
        ),
        list(
          "key" = "MISC_status-Convalescent", "name" = "Convalescent",
          "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
        )
      )
    ),
    "ScType-Spleen-human" = list(
      "key" = "67466668-bd95-11ed-9732-0242ac130003",
      "name" = "ScType-Spleen-human",
      "rootNode" = TRUE,
      "type" = "cellSets",
      "children" = list(
      list(
        "key" = "67466668-bd95-11ed-9732-0242ac130003", "name" = "Naive CD4+ T cells", "color" = "#77aadd",
        "cellIds" = unname(data$cells_id[1:(length(data$cells_id) / 2)])
      ),
      list(
        "key" = "67466668-bd95-11ed-9732-0242ac130003", "name" = "Platelets",
        "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
      )
    )
  )
)

  return(cell_sets)
}


with_fake_http(test_that("sendCellSetToApi sends patch request", {
  scr_cellset_object <- mock_scratchpad_cellset_object(5)
  expect_PATCH(sendCellsetToApi(scr_cellset_object, "api_url", "experiment_id", "cell_set_key", "auth"))
}))


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

  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)

  expect_true(data.table::is.data.table(parsed_cellsets))
  expect_equal(colnames(parsed_cellsets), expected_columns)
  expect_equal(length(children_cell_sets$louvain), length(unique(parsed_cellsets[cellset_type == "louvain", name])))
  expect_equal(sapply(children_cell_sets$louvain, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "louvain", name]))
  expect_equal(length(children_cell_sets$scratchpad), length(unique(parsed_cellsets[cellset_type == "scratchpad", name])))
  expect_equal(sapply(children_cell_sets$scratchpad, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "scratchpad", name]))
  expect_equal(length(children_cell_sets$sample), length(unique(parsed_cellsets[cellset_type == "sample", name])))
  expect_equal(sapply(children_cell_sets$sample, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "sample", name]))
  expect_equal(length(children_cell_sets$MISC_status), length(unique(parsed_cellsets[cellset_type == "MISC_status", name])))
  expect_equal(sapply(children_cell_sets$MISC_status, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "MISC_status", name]))
  expect_equal(length(children_cell_sets$`ScType-Spleen-human`), length(unique(parsed_cellsets[cellset_type == "ScType-Spleen-human", name])))
  expect_equal(sapply(children_cell_sets$`ScType-Spleen-human`, function(x) {
    x$name
  }), unique(parsed_cellsets[cellset_type == "ScType-Spleen-human", name]))
})

test_that("add_clusters adds cluster information as metadata columns to the seurat object", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)

  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  sctype_clusters <- parsed_cellsets[
    sapply(parsed_cellsets$cellset_type, function(cellset_type) {
      cell_sets[[cellset_type]]$type == "cellSets" &&
        cell_sets[[cellset_type]]$key != "louvain" &&
        cell_sets[[cellset_type]]$key != "scratchpad"
    }),
  ]

  sctype_clusters_list <- split(sctype_clusters, sctype_clusters[["cellset_type"]])
  for (sctype_group in sctype_clusters_list) {
    sctype_colname <- unique(sapply(sctype_group$cellset_type, function(x) {
      cell_sets[[x]]$name
    }))
    expect_true(sctype_colname %in% colnames(data@meta.data))
    sctype_dt <- sctype_group[, c("name", "cell_id")]
    data.table::setnames(sctype_dt, c(sctype_colname, "cells_id"))
    expect_equal(unique(data@meta.data[,sctype_colname]), unique(sctype_dt[, get(sctype_colname)]))
    expect_true(all.equal(sctype_dt, data@meta.data[,c(sctype_colname, "cells_id")],check.attributes = FALSE))
  }
})

test_that("add_clusters uses conventions that support re-upload with technology Seurat", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)

  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  # 'samples' column used for sample identity
  expect_true("samples" %in% colnames(data@meta.data))

  # active ident used as default louvain clusters
  expect_identical(as.character(data$seurat_clusters), as.character(Seurat::Idents(data)))
})


test_that("format_sctype_cell_sets correctly format cellset to be sent to the API", {
  data <- mock_scdata()
  cell_sets <- mock_cellset_from_python(data)
  species <- "human"
  tissue <- "Pancreas"
  active_assay <- "RNA"

  expected_cell_class_names <- c("key", "name", "rootNode", "type", "children")

  scale_data <- get_formatted_data(data, active_assay)
  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  data@meta.data <- suppressWarnings(
    run_sctype(scale_data, data@meta.data, tissue)
  )

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


test_that("get_feature_types correctly determines the ids_sym features type in the annot data frame", {
  data <- mock_scdata()
  annot <- data.table::as.data.table(data@misc$gene_annotations)
  annot <- annot[, .(input, original_name)]

  feature_types <- get_feature_types(annot)

  expect_equal(feature_types, IDS_SYM)
})


test_that("get_feature_types correctly determines the sym_ids features type in the annot data frame", {
  data <- mock_scdata()
  annot <- data.table::as.data.table(data@misc$gene_annotations)
  annot <- annot[, .(input, original_name)]

  annot_switched <- annot[,c(2,1)]

  feature_types <- get_feature_types(annot_switched)

  expect_equal(feature_types, SYM_IDS)
})


test_that("get_feature_types correctly determines the sym_sym features type in the annot data frame", {
  data <- mock_scdata()
  annot <- data.table::as.data.table(data@misc$gene_annotations)
  annot <- annot[, .(input, original_name)]

  annot_rep <- annot[,c(2,2)]

  feature_types <- get_feature_types(annot_rep)

  expect_equal(feature_types, SYM_SYM)
})


test_that("get_feature_types correctly determines the sym_sym features type in the annot data frame", {
  data <- mock_scdata()
  annot <- data.table::as.data.table(data@misc$gene_annotations)
  annot <- annot[, .(input, original_name)]

  annot_rep <- annot[,c(1,1)]

  feature_types <- get_feature_types(annot_rep)

  expect_equal(feature_types, IDS_IDS)
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


with_fake_http(test_that(
  "updateCellSetsThroughApi appends cellset when append = TRUE",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(10)
    req <- updateCellSetsThroughApi(
      scr_cellset_object,
      "api_url",
      "experiment_id",
      "cell_set_key",
      "auth",
      append = TRUE
    )

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    expect_equal(names(req_body[[2]]), "$append")
  }
))



with_fake_http(test_that(
  "updateCellSetsThroughApi prepend cellset when append = FALSE",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(10)
    req <- updateCellSetsThroughApi(
      scr_cellset_object,
      "api_url",
      "experiment_id",
      "cell_set_key",
      "auth",
      append = FALSE
    )

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    expect_equal(names(req_body[[2]]), "$prepend")
  }
))

test_that("get_quantile_cap (vectorized) produces same results as getQuantileCapOld for each column", {
  # Define the original non-vectorized implementation for comparison
  getQuantilCapOld <- function(x, quantile_threshold) {
    lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
    i <- 0.01
    while (lim == 0 && i + quantile_threshold <= 1) {
      lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
      i <- i + 0.01
    }
    return(lim)
  }
  
  # Create a test matrix with varied data patterns
  set.seed(42)
  
  # Create a sparse matrix with different characteristics per column
  test_matrix <- Matrix::Matrix(
    c(
      0, 1, 2, 3, 4, 5,           # Column 1: typical data
      0, 0, 0, 1, 2, 10,          # Column 2: mostly zeros then spikes
      5, 6, 7, 8, 9, 10,          # Column 3: uniform data
      0, 0, 0, 0, 0, 100,         # Column 4: zeros with one large value
      1, 2, 3, 4, 5, 6            # Column 5: linear progression
    ),
    nrow = 6,
    ncol = 5,
    sparse = TRUE
  )
  
  # Test with different quantile thresholds
  quantile_thresholds <- c(0.75, 0.9, 0.95)
  
  for (threshold in quantile_thresholds) {
    # Get vectorized results (all columns at once)
    vectorized_caps <- get_quantile_cap(test_matrix, threshold)
    
    # Get non-vectorized results (column by column)
    old_caps <- sapply(
      1:ncol(test_matrix),
      function(i) getQuantilCapOld(test_matrix[, i], threshold)
    )
    
    # Compare results
    expect_equal(
      vectorized_caps,
      as.numeric(old_caps),
      tolerance = 1e-10,
      label = sprintf("Quantile threshold: %s", threshold)
    )
  }
})

test_that("get_quantile_cap (vectorized) handles edge cases like old implementation", {
  # Define the original non-vectorized implementation for comparison
  getQuantilCapOld <- function(x, quantile_threshold) {
    lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
    i <- 0.01
    while (lim == 0 && i + quantile_threshold <= 1) {
      lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
      i <- i + 0.01
    }
    return(lim)
  }
  
  # Test with all zeros
  zero_matrix <- Matrix::Matrix(rep(0, 12), nrow = 6, ncol = 2, sparse = TRUE)
  
  vectorized_result <- get_quantile_cap(zero_matrix, 0.75)
  old_result <- c(
    getQuantilCapOld(zero_matrix[, 1], 0.75),
    getQuantilCapOld(zero_matrix[, 2], 0.75)
  )
  
  expect_equal(vectorized_result, as.numeric(old_result))
  
  # Test with NAs
  na_matrix <- Matrix::Matrix(
    c(NA, 1, 2, NA, 3, 4, NA, 5, 6),
    nrow = 3,
    ncol = 3,
    sparse = TRUE
  )
  
  vectorized_result <- get_quantile_cap(na_matrix, 0.75)
  old_result <- c(
    getQuantilCapOld(na_matrix[, 1], 0.75),
    getQuantilCapOld(na_matrix[, 2], 0.75),
    getQuantilCapOld(na_matrix[, 3], 0.75)
  )
  
  expect_equal(vectorized_result, as.numeric(old_result))
})

test_that("get_quantile_cap increments threshold when initial quantile is zero", {
  # Create a matrix where high quantiles are zero but low values exist
  # High quantile thresholds on sparse data can return 0
  test_matrix <- Matrix::Matrix(
    c(
      0, 0, 0, 0, 0, 1,        # Column 1: mostly zeros with one value at 1
      0, 0, 0, 0, 0, 0.5,      # Column 2: mostly zeros with one value at 0.5
      0, 0, 0, 0, 0, 2         # Column 3: mostly zeros with one value at 2
    ),
    nrow = 6,
    ncol = 3,
    sparse = TRUE
  )
  
  # At 0.95 quantile, all columns would return 0 initially
  result <- get_quantile_cap(test_matrix, 0.95)
  
  # Should not return 0 for columns with non-zero values
  expect_true(all(result != 0))
  
  # The result should be close to the actual non-zero values
  expect_true(result[1] > 0 && result[1] <= 1)
  expect_true(result[2] > 0 && result[2] <= 0.5)
  expect_true(result[3] > 0 && result[3] <= 2)
})

test_that("get_quantile_cap returns 0 when all values are zero", {
  # Create a matrix with all zeros
  zero_matrix <- Matrix::Matrix(rep(0, 20), nrow = 5, ncol = 4, sparse = TRUE)
  
  result <- get_quantile_cap(zero_matrix, 0.75)
  
  # Should return 0 for all columns when all values are zero
  expect_equal(result, rep(0, ncol(zero_matrix)))
})

test_that("get_quantile_cap handles high quantiles correctly", {
  # Create a test matrix with a range of values
  test_matrix <- Matrix::Matrix(
    c(
      1, 2, 3, 4, 5, 6,            # Column 1: 1-6
      0.1, 0.2, 0.3, 0.4, 0.5, 10  # Column 2: mostly small values with one large
    ),
    nrow = 6,
    ncol = 2,
    sparse = TRUE
  )
  
  # Get quantiles at different thresholds
  result_75 <- get_quantile_cap(test_matrix, 0.75)
  result_90 <- get_quantile_cap(test_matrix, 0.90)
  result_99 <- get_quantile_cap(test_matrix, 0.99)
  
  # Higher quantiles should generally give higher or equal values
  expect_true(result_90[1] >= result_75[1])
  expect_true(result_99[1] >= result_90[1])
  
  # For column with sparse values, even high quantiles should not return 0 if any non-zero exists
  expect_true(result_99[2] > 0)
})
