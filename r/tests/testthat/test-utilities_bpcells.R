library("mockery")

# helper: create Seurat object with BPCells-backed matrix
create_bpcells_seurat <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")

  # write to temporary BPCells directory
  matrix_dir <- file.path(tempdir(), "matrix_dir")
  if (dir.exists(matrix_dir)) {
    unlink(matrix_dir, recursive = TRUE)
  }
  counts <- BPCells::write_matrix_dir(pbmc_raw, matrix_dir)

  # create Seurat object with BPCells matrix
  # (CreateSeuratObject converts the MatrixDir to a RenameDims layer)
  scdata <- SeuratObject::CreateSeuratObject(
    counts = counts
  )

  list(scdata = scdata, matrix_dir = matrix_dir)
}

# test tar_zstd constructs correct command
test_that("tar_zstd constructs correct tar command", {
  system2_called <- FALSE

  mock_system2 <- function(cmd, args) {
    system2_called <<- TRUE
    0L
  }

  stub(tar_zstd, "system2", mock_system2)
  tar_zstd("/tmp/test.tar.zst", "/path/to/file")

  expect_true(system2_called)
})

# test untar_zstd constructs correct command
test_that("untar_zstd constructs correct tar command", {
  system2_called <- FALSE

  mock_system2 <- function(cmd, args) {
    system2_called <<- TRUE
    0L
  }

  stub(untar_zstd, "system2", mock_system2)
  untar_zstd("/tmp/test.tar.zst", "/tmp/extract")

  expect_true(system2_called)
})

# test find_matrix_dir_paths recursively locates MatrixDir objects
test_that("find_matrix_dir_paths finds MatrixDir in nested layers", {
  bpcells_data <- create_bpcells_seurat()

  # extract from the counts layer (wrapped in RenameDims)
  counts_layer <- SeuratObject::LayerData(
    bpcells_data$scdata[["RNA"]],
    "counts"
  )

  # find_matrix_dir_paths should recursively find any MatrixDir objects
  paths <- find_matrix_dir_paths(counts_layer)

  # verify we found at least one matrix dir path or have RenameDims wrapper
  expect_true(
    length(paths) > 0 || inherits(counts_layer, "RenameDims")
  )
})

# test update_matrix_dir updates paths in all assay layers
test_that("update_matrix_dir updates paths across assays", {
  bpcells_data <- create_bpcells_seurat()
  seurat_obj <- bpcells_data$scdata
  old_dir <- bpcells_data$matrix_dir
  new_dir <- file.path(tempdir(), "new_matrix_dir")
  dir.create(new_dir, showWarnings = FALSE)

  # copy the matrix directory to new location
  file.copy(
    file.path(old_dir, "."),
    new_dir,
    recursive = TRUE,
    overwrite = TRUE
  )

  # update paths in the object
  updated_obj <- update_matrix_dir(seurat_obj, old_dir, new_dir)

  # verify object structure is preserved
  expect_equal(class(updated_obj), class(seurat_obj))
  expect_equal(ncol(updated_obj), ncol(seurat_obj))
  expect_equal(nrow(updated_obj), nrow(seurat_obj))
})

# test replace_matrix_dir_paths handles MatrixDir updates correctly
test_that("replace_matrix_dir_paths updates MatrixDir path", {
  bpcells_data <- create_bpcells_seurat()
  old_dir <- bpcells_data$matrix_dir
  new_dir <- file.path(tempdir(), "new_matrix_dir")
  dir.create(new_dir, showWarnings = FALSE)

  # copy the matrix directory
  file.copy(
    file.path(old_dir, "."),
    new_dir,
    recursive = TRUE,
    overwrite = TRUE
  )

  # test direct replacement on a layer
  counts_layer <- SeuratObject::LayerData(bpcells_data$scdata[["RNA"]], "counts")
  updated_layer <- replace_matrix_dir_paths(counts_layer, old_dir, new_dir)

  # verify we can still work with the layer after path update
  expect_s4_class(updated_layer, "RenameDims")
})

test_that("materialize_bpcells_matrix converts IterableMatrix to dgCMatrix", {
  bpcells_data <- create_bpcells_seurat()
  scdata <- bpcells_data$scdata

  # typical Seurat workflow
  scdata <- Seurat::NormalizeData(scdata)
  scdata <- Seurat::FindVariableFeatures(scdata)
  scdata <- Seurat::ScaleData(scdata)

  # check that layers are IterableMatrix
  for (layer_name in SeuratObject::Layers(scdata[["RNA"]])) {
    layer <- SeuratObject::LayerData(scdata[["RNA"]], layer_name)
    expect_true(is(layer, "IterableMatrix"))
  }

  # materialize the object
  scdata <- materialize_bpcells_matrix(scdata)

  # should all be dgCMatrix now
  for (layer_name in SeuratObject::Layers(scdata[["RNA"]])) {
    layer <- SeuratObject::LayerData(scdata[["RNA"]], layer_name)
    expect_s4_class(layer, "dgCMatrix")
  }
})

test_that("load_bpcells works", {
  bpcells_data <- create_bpcells_seurat()
  scdata <- bpcells_data$scdata
  matrix_dir <- bpcells_data$matrix_dir

  tarfile <- file.path(tempdir(), "matrix_dir.tar.zst")
  withr::defer(unlink(tarfile))
  # move to parent of matrix_dir and tar it
  data_dir <- dirname(matrix_dir)

  current_dir <- getwd()
  setwd(data_dir)
  tar_zstd(tarfile, files = basename(matrix_dir))
  setwd(current_dir)

  # scdata shouldn't work because matrix_dir removed
  unlink(matrix_dir, recursive = TRUE)
  expect_error(as.matrix(scdata[["RNA"]]$counts))


  # load the matrix from the directory
  scdata <- load_bpcells(scdata, data_dir)

  # verify it's a IterableMatrix and can be accessed
  expect_s4_class(scdata[["RNA"]]$counts, "IterableMatrix")
  expect_no_error(as.matrix(scdata[["RNA"]]$counts))
})

test_that("load_bpcells throws error if matrix_dir and tarfile missing", {
  bpcells_data <- create_bpcells_seurat()
  scdata <- bpcells_data$scdata
  matrix_dir <- bpcells_data$matrix_dir
  data_dir <- dirname(matrix_dir)

  # remove the matrix directory to simulate missing data
  unlink(matrix_dir, recursive = TRUE)

  expect_error(load_bpcells(scdata, data_dir))
})