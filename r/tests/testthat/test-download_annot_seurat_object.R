
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

# write the download to a temp zip instead of DOWNLOAD_ZIP_PATH
download_zip <- file.path(tempdir(), "processed.zip")

mockery::stub(
  DownloadAnnotSeuratObject,
  "save_download_zip",
  function(data, zip_path) save_download_zip(data, download_zip)
)

stubbed_DownloadAnnotSeuratObject <- function(req, data) {
  DownloadAnnotSeuratObject(req, data)
}

# extract the download zip to a new folder, as a user would
unzip_download <- function() {
  out_dir <- withr::local_tempdir(.local_envir = parent.frame())
  utils::unzip(download_zip, exdir = out_dir)
  out_dir
}

# run the zip's load_object.R from the folder it was extracted to
run_load_script <- function(out_dir) {
  env <- new.env()
  withr::with_dir(out_dir, sys.source("load_object.R", envir = env))
  env$scdata
}

test_that("DownloadAnnotSeuratObject saves a zip with the Seurat object and a load script", {
  data <- mock_scdata(with_umap = TRUE)
  req <- mock_req(data)

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))
  expect_equal(res, DOWNLOAD_ZIP_PATH)

  out_dir <- unzip_download()
  # in-memory counts: no matrix dir to ship
  expect_setequal(list.files(out_dir), c("object.rds", "load_object.R"))

  loaded_data <- readRDS(file.path(out_dir, "object.rds"))
  expect_true(is(loaded_data, "Seurat"))
  expect_true("seurat_clusters" %in% colnames(loaded_data@meta.data))
  expect_true("umap" %in% names(loaded_data@reductions))
})


test_that("DownloadAnnotSeuratObject works with Seurat projects", {
  data <- mock_scdata(with_umap = TRUE)
  req <- mock_req(data)
  req$isObj2s <- TRUE
  req$embedding <- NULL

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_equal(res, DOWNLOAD_ZIP_PATH)
  expect_true(file.exists(download_zip))
})


test_that("DownloadAnnotSeuratObject keeps BPCells counts disk-backed with their matrix dir", {
  data <- mock_scdata(with_umap = TRUE, use_bpcells = TRUE)
  req <- mock_req(data)

  suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  out_dir <- unzip_download()
  expect_setequal(list.files(out_dir), c("object.rds", "load_object.R", "matrix_dir"))
  expect_setequal(list.files(file.path(out_dir, "matrix_dir")), list.files(get_matrix_dirs(data)))

  # counts are not read into memory
  loaded_data <- readRDS(file.path(out_dir, "object.rds"))
  expect_true(is(loaded_data[["RNA"]]$counts, "IterableMatrix"))
})


test_that("the download load script loads the disk-backed object from where it was extracted", {
  data <- mock_scdata(with_umap = TRUE, use_bpcells = TRUE)
  req <- mock_req(data)
  expected_counts <- as(data[["RNA"]]$counts, "dgCMatrix")

  suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))
  out_dir <- unzip_download()

  # the worker's matrix dir is gone: the counts must come from the zip's matrix_dir
  unlink(get_matrix_dirs(data), recursive = TRUE)

  scdata <- run_load_script(out_dir)

  expect_true(is(scdata, "Seurat"))
  expect_equal(
    normalizePath(get_matrix_dirs(scdata)),
    normalizePath(file.path(out_dir, "matrix_dir"))
  )
  expect_equal(as(scdata[["RNA"]]$counts, "dgCMatrix"), expected_counts)
})
