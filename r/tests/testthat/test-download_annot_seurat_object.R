mock_color_pool <- function(n) {
  paste0("color_", 1:n)
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
  pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")
  pbmc_small <- SeuratObject::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  pbmc_small@misc$color_pool <- mock_color_pool(20)

  pbmc_small <- Seurat::NormalizeData(
    pbmc_small,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)

  # scale and PCA
  pbmc_small <- Seurat::NormalizeData(
    pbmc_small,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::RunPCA(pbmc_small, verbose = FALSE, npcs = 10)
  pbmc_small@misc[["active.reduction"]] <- "pca"

  # run UMAP
  pbmc_small <- suppressWarnings(
    Seurat::RunUMAP(pbmc_small, dims = 1:10, verbose = FALSE)
  )

  # add sample metadata
  pbmc_small@meta.data$samples <- rep(
    "0000-0000-0000-0000-0000",
    nrow(pbmc_small@meta.data)
  )

  return(pbmc_small)
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
  data <- mock_scdata()
  req <- mock_req(data)

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)
})


test_that("DownloadAnnotSeuratObject works with Seurat projects", {
  data <- mock_scdata()
  req <- mock_req(data)
  req$isObj2s <- TRUE
  req$embedding <- NULL

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)
})

# create BPCells-backed Seurat object with same structure as mock_scdata
mock_bpcells_scdata <- function() {
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
  pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")

  # write to temporary BPCells directory
  temp_bpcells_dir <- file.path(tempdir(), "mock_bpcells_matrix")
  if (dir.exists(temp_bpcells_dir)) {
    unlink(temp_bpcells_dir, recursive = TRUE)
  }
  BPCells::write_matrix_dir(pbmc_raw, temp_bpcells_dir)

  # create Seurat object from BPCells matrix
  pbmc_small <- SeuratObject::CreateSeuratObject(
    counts = BPCells::open_matrix_dir(temp_bpcells_dir)
  )

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  pbmc_small@misc$color_pool <- mock_color_pool(20)

  pbmc_small <- Seurat::NormalizeData(
    pbmc_small,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)

  # scale and PCA
  pbmc_small <- Seurat::RunPCA(pbmc_small, verbose = FALSE, npcs = 10)
  pbmc_small@misc[["active.reduction"]] <- "pca"

  # run UMAP
  pbmc_small <- suppressWarnings(Seurat::RunUMAP(
    pbmc_small,
    dims = 1:10,
    verbose = FALSE
  ))

  # add sample metadata
  pbmc_small@meta.data$samples <- rep(
    "0000-0000-0000-0000-0000",
    nrow(pbmc_small@meta.data)
  )

  # store BPCells directory path for cleanup later
  pbmc_small@misc$bpcells_dir <- temp_bpcells_dir

  return(pbmc_small)
}

test_that("DownloadAnnotSeuratObject works with BPCells disk-backed objects", {
  data <- mock_bpcells_scdata()
  req <- mock_req(data)

  res <- suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, data))

  expect_type(res, "character")
  expect_equal(res, RDS_PATH)

  # cleanup
  bpcells_dir <- data@misc$bpcells_dir
  if (dir.exists(bpcells_dir)) {
    unlink(bpcells_dir, recursive = TRUE)
  }
})

test_that("DownloadAnnotSeuratObject maintains BPCells functionality after save/load", {
  data <- mock_bpcells_scdata()
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

  # cleanup
  bpcells_dir <- data@misc$bpcells_dir
  if (dir.exists(bpcells_dir)) {
    unlink(bpcells_dir, recursive = TRUE)
  }
})

test_that("DownloadAnnotSeuratObject works with tar.zst BPCells workflow", {
  # simulate the real worker workflow where BPCells comes in tar.zst format
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")

  # create temporary directories simulating data flow
  data_dir_1 <- file.path(tempdir(), "data_extraction_dir")
  data_dir_2 <- file.path(tempdir(), "data_result_dir")
  on.exit({
    if (dir.exists(data_dir_1)) unlink(data_dir_1, recursive = TRUE)
    if (dir.exists(data_dir_2)) unlink(data_dir_2, recursive = TRUE)
  })

  if (!dir.exists(data_dir_1)) dir.create(data_dir_1, recursive = TRUE)
  if (!dir.exists(data_dir_2)) dir.create(data_dir_2, recursive = TRUE)

  # write matrix to tar.zst in first directory
  matrix_dir_1 <- file.path(data_dir_1, "matrix_dir")
  BPCells::write_matrix_dir(pbmc_raw, matrix_dir_1)

  tarfile <- file.path(data_dir_1, "matrix_dir.tar.zst")
  system(paste("tar --zstd -cf", tarfile, "-C", data_dir_1, "matrix_dir"))

  # simulate loading: extract tar, create Seurat object
  untar_zstd(tarfile, exdir = data_dir_1)
  matrix_dir_path_1 <- file.path(data_dir_1, "matrix_dir")

  scdata <- SeuratObject::CreateSeuratObject(
    counts = BPCells::open_matrix_dir(matrix_dir_path_1)
  )
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # normalize and add UMAP
  scdata <- Seurat::NormalizeData(scdata, verbose = FALSE)
  scdata <- Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  scdata <- Seurat::RunPCA(scdata, verbose = FALSE, npcs = 10)
  scdata <- suppressWarnings(Seurat::RunUMAP(scdata, dims = 1:10, verbose = FALSE))

  # prepare request and embedding
  req <- list(
    body = list(
      cellSets = mock_cellset_from_python(scdata),
      embedding = mock_embedding_data(scdata),
      embeddingMethod = "umap",
      isObj2s = FALSE
    )
  )

  # call DownloadAnnotSeuratObject
  res <- tryCatch(
    {
      suppressWarnings(stubbed_DownloadAnnotSeuratObject(req, scdata))
    },
    error = function(e) {
      stop("DownloadAnnotSeuratObject failed with BPCells: ", e$message)
    }
  )

  expect_type(res, "character")

  # verify file was saved
  expect_true(file.exists(file.path(tempdir(), "data", "processed.rds")))
})

test_that(
  "DownloadAnnotSeuratObject converts IterableMatrix to dgCMatrix",
  {
    # create BPCells object
    data <- mock_bpcells_scdata()
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

    # get the BPCells directory path before cleanup
    bpcells_dir <- data@misc$bpcells_dir

    # load the saved file
    loaded_data <- readRDS(file.path(tempdir(), "data", "processed.rds"))

    # verify loaded object exists
    expect_true(is(loaded_data, "Seurat"))

    # now delete the BPCells directory to ensure matrices aren't disk-backed
    unlink(bpcells_dir, recursive = TRUE)

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
    expect_true(nrow(loaded_data@assays$RNA$counts) > 0)
    expect_true(ncol(loaded_data@assays$RNA$counts) > 0)
    expect_true(nrow(loaded_data@assays$RNA$data) > 0)
    expect_true(ncol(loaded_data@assays$RNA$data) > 0)
    expect_true(nrow(loaded_data@assays$RNA$scale.data) > 0)
    expect_true(ncol(loaded_data@assays$RNA$scale.data) > 0)

    # cleanup
    if (dir.exists(bpcells_dir)) {
      unlink(bpcells_dir, recursive = TRUE)
    }
  }
)
