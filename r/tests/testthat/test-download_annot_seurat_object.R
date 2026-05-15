
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
          "cellIds" = unname(data$cells_id[1:ceiling(length(data$cells_id)/3)])
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


mock_embedding_data <- function(data) {
  embedding_data <- data@reductions$umap@cell.embeddings
  embedding_data <- unname(embedding_data)
  embedding_data <- lapply(
    seq_len(nrow(embedding_data)),
    function(i) embedding_data[i, ]
  )
  return(embedding_data)
}

mock_req <- function(data, reduction = "umap") {
  req <- list(
    body = list(
      cellSets = mock_cellset_from_python(data),
      embedding = mock_embedding_data(data),
      embeddingMethod = reduction,
      isObj2s = FALSE
    )
  )

  return(req)
}

stub_saveRDS <- function(data, fpath) {
  rdir <- file.path(tempdir(), "data")
  if (!dir.exists(rdir)) {
    dir.create(rdir, recursive = TRUE)
  }
  fpath <- file.path(rdir, "processed.rds")
  saveRDS(data, fpath)
}

mockery::stub(
  DownloadAnnotSeuratObject,
  "saveRDS",
  stub_saveRDS
)

stubbed_DownloadAnnotSeuratObject <- function(req, data) {
  DownloadAnnotSeuratObject(req, data)
}

test_that("DownloadAnnotSeuratObject saves the Seurat object using the correct path", {
  data <- mock_scdata(with_umap = TRUE)
  req <- mock_req(data)

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)
})


test_that("DownloadAnnotSeuratObject works with Seurat projects", {
  data <- mock_scdata(with_umap = TRUE)
  req <- mock_req(data)
  req$isObj2s <- TRUE
  req$embedding <- NULL

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)
})


test_that("DownloadAnnotSeuratObject works with BPCells disk-backed objects", {
  data <- mock_scdata(with_umap = TRUE, use_bpcells = TRUE)
  req <- mock_req(data)

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)
})

test_that("DownloadAnnotSeuratObject maintains BPCells functionality after save/load", {
  data <- mock_scdata(with_umap = TRUE, use_bpcells = TRUE)
  req <- mock_req(data)

  # call DownloadAnnotSeuratObject which will save the file
  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  # attempt to load the saved RDS file
  loaded_data <- tryCatch(
    {
      readRDS(file.path(tempdir(), "data", "processed.rds"))
    },
    error = function(e) {
      NULL
    }
  )

  # verify the loaded data maintains structure
  expect_false(is.null(loaded_data))
  expect_true(is(loaded_data, "Seurat"))

  # check that we can access metadata that was added by add_clusters
  expect_true("seurat_clusters" %in% colnames(loaded_data@meta.data))

  # check that we can access the UMAP embedding that was added
  expect_true("umap" %in% names(loaded_data@reductions))
})


test_that("DownloadAnnotSeuratObject converts IterableMatrix to dgCMatrix", {
    # create BPCells object
    data <- mock_scdata(with_umap = TRUE, use_bpcells = TRUE)
    req <- mock_req(data)

    # verify the original object has IterableMatrix
    expect_true(
      is(data@assays$RNA$counts, "IterableMatrix"),
      label = "Original should have IterableMatrix"
    )

    # call DownloadAnnotSeuratObject which should save converted data
    res <- suppressWarnings(
      stubbed_DownloadAnnotSeuratObject(req, data)
    )

    # load the saved file
    loaded_data <- readRDS(file.path(tempdir(), "data", "processed.rds"))

    # verify loaded object exists
    expect_true(is(loaded_data, "Seurat"))

    # delete the original matrix_dir to ensure matrices aren't disk-backed
    matrix_dir <- get_matrix_dirs(data)

    expect_true(dir.exists(matrix_dir))
    unlink(matrix_dir, recursive = TRUE)
    expect_false(dir.exists(matrix_dir))

    # verify all matrix slots are accessible and are dgCMatrix
    # (not IterableMatrix which would require disk access)

    # Check counts matrix
    expect_true(
      is(loaded_data@assays$RNA$counts, "dgCMatrix"),
      label = "counts should be dgCMatrix after load"
    )

    # Check data matrix (normalized)
    expect_true(
      is(loaded_data@assays$RNA$data, "dgCMatrix"),
      label = "data should be dgCMatrix after load"
    )

    # Check scale.data matrix
    expect_true(
      is(loaded_data@assays$RNA$scale.data, "dgCMatrix"),
      label = "scale.data should be dgCMatrix after load"
    )

    # verify we can actually access data without errors
    expect_no_error(loaded_data@assays$RNA$counts[1:10, 1:10])
    expect_no_error(loaded_data@assays$RNA$data[1:10, 1:10])
    expect_no_error(loaded_data@assays$RNA$scale.data[1:10, 1:10])
  }
)
