mock_cellset_from_python <- function(data) {
  cell_sets <- list(
    "louvain" = list(
      "key" = "louvain",
      "name" = "louvain clusters",
      "rootNode" = TRUE,
      "type" = "cellSets",
      "children" = list(
        list(
          "key" = "louvain-0",
          "name" = "Cluster 0",
          "rootNode" = FALSE,
          "type" = "cellSets",
          "color" = "#77aadd",
          "cellIds" = unname(
            data$cells_id[1:ceiling(length(data$cells_id) / 3)]
          )
        ),
        list(
          "key" = "louvain-1",
          "name" = "Cluster 1",
          "rootNode" = FALSE,
          "type" = "cellSets",
          "color" = "#77aadd",
          "cellIds" = unname(
            data$cells_id[
              (ceiling(length(data$cells_id) / 3 + 1)):
                (ceiling(length(data$cells_id) / 3 * 2))
            ]
          )
        ),
        list(
          "key" = "louvain-2",
          "name" = "Cluster 2",
          "rootNode" = FALSE,
          "type" = "cellSets",
          "color" = "#ee8866",
          "cellIds" = unname(
            data$cells_id[
              (ceiling(length(data$cells_id) / 3 * 2 + 1)):
                (length(data$cells_id))
              ]
            )
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
          "key" = "scratchpad-0",
          "name" = "Custom 0",
          "rootNode" = FALSE,
          "type" = "cellSets",
          "color" = "#77aadd",
          "cellIds" = unname(sample(data$cells_id, 5))
        ),
        list(
          "key" = "scratchpad-1",
          "name" = "Custom 1",
          "rootNode" = FALSE,
          "type" = "cellSets",
          "color" = "#ee8866",
          "cellIds" = unname(sample(data$cells_id, 10))
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
          "key" = "a636ec18-4ba3-475b-989d-0a5b2",
          "name" = "P13 Acute MISC",
          "color" = "#77aadd",
          "cellIds" = unname(data$cells_id[1:(length(data$cells_id) / 2)])
        ),
        list(
          "key" = "eec701f2-5762-4b4f-953d-6aba8",
          "name" = "P13 Convalescent MISC",
          "color" = "#ee8866",
          "cellIds" = unname(
            data$cells_id[
              (length(data$cells_id) / 2 + 1):(length(data$cells_id))
            ]
          )
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
          "key" = "MISC_status-Acute",
          "name" = "Acute",
          "color" = "#77aadd",
          "cellIds" = unname(data$cells_id[1:(length(data$cells_id) / 2)])
        ),
        list(
          "key" = "MISC_status-Convalescent",
          "name" = "Convalescent",
          "color" = "#ee8866",
          "cellIds" = unname(
            data$cells_id[
              (length(data$cells_id) / 2 + 1):(length(data$cells_id))
            ]
          )
        )
      )
    )
  )

  return(cell_sets)
}

mock_req <- function(data) {
  req <- list(
    body = list(
      cellSets = mock_cellset_from_python(data),
      species = "human",
      tissue = "Pancreas"
    )
  )

  return(req)
}


stub_updateCellSetsThroughApi <- function(
  cell_sets_object,
  api_url,
  experiment_id,
  cell_set_key,
  auth_JWT,
  append = TRUE
) {

  # empty function to simplify mocking. we test patching independently.
}


stubbed_ScTypeAnnotate <- function(req, data) {
  mockery::stub(ScTypeAnnotate,
                "updateCellSetsThroughApi",
                stub_updateCellSetsThroughApi)
  ScTypeAnnotate(req, data)
}


test_that("add_gene_symbols adds gene symbols to the count matrix", {
  data <- mock_scdata()
  active_assay <- "RNA"

  scale_data <- data.table::as.data.table(
    data[[active_assay]]$scale.data,
    keep.rownames = "input"
  )
  scale_data <- add_gene_symbols(scale_data, data)

  expect_equal(colnames(scale_data)[2], "original_name")
})


test_that("add_gene_symbols produces an error if there are no gene symbols in the annot data frame", {
  data <- mock_scdata()
  active_assay <- "RNA"

  scale_data <- data.table::as.data.table(
    data[[active_assay]]$scale.data,
    keep.rownames = "input"
  )

  annot <- data.table::as.data.table(data@misc$gene_annotations)

  annot_mod <- annot
  annot_mod$name <- annot_mod$input
  annot_mod$original_name <- annot_mod$input
  data@misc$gene_annotations <- annot_mod

  expect_error(
    add_gene_symbols(scale_data, data),
    "Features file doesn't contain gene symbols."
  )
})


test_that("collapse_genes collapses duplicated gene symbols", {
  data <- mock_scdata()
  active_assay <- "RNA"

  scale_data <- data.table::as.data.table(
    data[[active_assay]]$scale.data,
    keep.rownames = "input"
  )
  scale_data <- add_gene_symbols(scale_data, data)

  # duplicate first 5 rows
  scale_data_dup <- rbind(scale_data, scale_data[1:5, ])

  scale_data_dedup <- collapse_genes(scale_data_dup)

  expect_true(nrow(scale_data_dup) != nrow(scale_data_dedup))
  expect_true(!any(duplicated(scale_data_dedup$original_name)))
})


test_that("format_matrix produces a matrix in the expected format", {
  data <- mock_scdata()
  active_assay <- "RNA"

  scale_data <- data.table::as.data.table(
    data[[active_assay]]$scale.data,
    keep.rownames = "input"
  )
  scale_data <- add_gene_symbols(scale_data, data)
  scale_data <- collapse_genes(scale_data)

  scale_data_formatted <- format_matrix(scale_data)

  expect_type(scale_data_formatted, "double")
  expect_true(!"input" %in% colnames(scale_data_formatted))
  expect_identical(rownames(scale_data_formatted), scale_data$original_name)
})


test_that("run_sctype adds cluster annotations as a new metadata column", {
  data <- mock_scdata()
  active_assay <- "RNA"
  req <- mock_req(data)
  tissue <- req$body$tissue
  species <- req$body$species
  cell_sets <- req$body$cellSets

  scale_data <- get_formatted_data(data, active_assay)
  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  data@meta.data <- suppressWarnings(
    run_sctype(scale_data, data@meta.data, tissue)
  )

  expect_true("customclassif" %in% colnames(data@meta.data))
  expect_equal(all(is.na(data@meta.data$customclassif)), FALSE)
})


test_that("run_sctype produces correct snapshots", {
  data <- mock_scdata()
  active_assay <- "RNA"
  req <- mock_req(data)
  tissue <- req$body$tissue
  species <- req$body$species
  cell_sets <- req$body$cellSets

  scale_data <- get_formatted_data(data, active_assay)
  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  data@meta.data <- suppressWarnings(
    run_sctype(scale_data, data@meta.data, tissue)
  )

  expect_snapshot(data@meta.data[c("cells_id","customclassif")])
})


test_that("ScTypeAnnotate produces correct annotations", {
  data <- mock_scdata()
  active_assay <- "RNA"
  req <- mock_req(data)

  # fix sampled colors
  set.seed(1)
  res <- suppressWarnings(stubbed_ScTypeAnnotate(req, data))

  # clean up time-based stuff
  mock_uuid <- "this_is_not_an_uuid"
  res$key <- mock_uuid
  for (i in 1:length(res$children)) {
    res$children[[i]]$key <- mock_uuid
  }

  expect_snapshot(res)
})

test_that("ScTypeAnnotate works with bpcells", {
  data <- mock_scdata(use_bpcells = TRUE)
  active_assay <- "RNA"
  req <- mock_req(data)

  # fix sampled colors
  expect_no_error(suppressWarnings(
    stubbed_ScTypeAnnotate(req, data)
  ))
})


test_that("ScTypeAnnotate works with sketch assay", {
  data <- mock_scdata(use_bpcells = TRUE, nreps = 10)
  data <- suppressWarnings(mock_sketch(data))
  active_assay <- "RNA"
  req <- mock_req(data)

  # fix sampled colors
  set.seed(1)
  res <- expect_no_error(suppressWarnings(
    stubbed_ScTypeAnnotate(req, data)
  ))
})
