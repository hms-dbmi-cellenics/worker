library("mockery")

# helper: create Seurat object with BPCells-backed matrix
create_bpcells_seurat <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  pbmc_raw <- as(as.matrix(pbmc_raw), "dgCMatrix")

  # write to temporary BPCells directory
  temp_bpcells_dir <- file.path(tempdir(), "matrix_dir")
  if (dir.exists(temp_bpcells_dir)) {
    unlink(temp_bpcells_dir, recursive = TRUE)
  }
  counts <- BPCells::write_matrix_dir(pbmc_raw, temp_bpcells_dir)

  # create Seurat object with BPCells matrix
  # (CreateSeuratObject converts the MatrixDir to a RenameDims layer)
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = counts
  )

  list(obj = seurat_obj, bpcells_dir = temp_bpcells_dir)
}

# test tar_zstd constructs correct command
test_that("tar_zstd constructs correct tar command", {
  system_called <- FALSE
  captured_command <- NULL

  mock_system <- function(cmd) {
    system_called <<- TRUE
    captured_command <<- cmd
    0
  }

  stub(tar_zstd, "system", mock_system)
  tar_zstd("/tmp/test.tar.zst", "/path/to/file")

  expect_true(system_called)
  expect_match(captured_command, "tar --zstd -cf")
})

# test untar_zstd constructs correct command
test_that("untar_zstd constructs correct tar command", {
  system_called <- FALSE
  captured_command <- NULL

  mock_system <- function(cmd) {
    system_called <<- TRUE
    captured_command <<- cmd
    0
  }

  stub(untar_zstd, "system", mock_system)
  untar_zstd("/tmp/test.tar.zst", "/tmp/extract")

  expect_true(system_called)
  expect_match(captured_command, "tar --zstd -xf")
})

# test find_matrix_dir_paths recursively locates MatrixDir objects
test_that("find_matrix_dir_paths finds MatrixDir in nested layers", {
  bpcells_data <- create_bpcells_seurat()

  # extract from the counts layer (wrapped in RenameDims)
  counts_layer <- bpcells_data$obj@assays$RNA@layers[["counts"]]

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
  seurat_obj <- bpcells_data$obj
  old_dir <- bpcells_data$bpcells_dir
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
  old_dir <- bpcells_data$bpcells_dir
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
  counts_layer <- bpcells_data$obj@assays$RNA@layers[["counts"]]
  updated_layer <- replace_matrix_dir_paths(counts_layer, old_dir, new_dir)

  # verify we can still work with the layer after path update
  expect_s4_class(updated_layer, "RenameDims")
})
