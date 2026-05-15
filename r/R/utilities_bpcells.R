load_bpcells <- function(data, data_dir) {
  # untar tarfiles
  tarfile <- list.files(
    data_dir,
    pattern = "matrix_dir[.]tar[.]zst",
    full.names = TRUE
  )

  if (!dir.exists(file.path(data_dir, "matrix_dir"))) {
    # need to throw error otherwise worker will report success
    # but the data won't be usable since the matrix_dir won't be extracted
    if (!length(tarfile)) stop("matrix_dir.tar.zst not found")

    message("Extracting BPCells matrix tarfile...")
    untar_zstd(tarfile, exdir = data_dir)
  }

  old_dir <- get_matrix_dirs(data)
  new_dir <- file.path(data_dir, basename(old_dir))
  data <- update_matrix_dir(data, old_dir, new_dir)

  return(data)
}


# TODO: newer R utils::tar supports zstd
tar_zstd <- function(tarfile, files) {
  status <- system2("tar", c("--zstd", "-cf", tarfile, files))
  if (!identical(status, 0L)) {
    stop("Failed to create zstd tar archive: ", tarfile)
  }
}

untar_zstd <- function(tarfile, exdir) {
  status <- system2("tar", c("--zstd", "-xf", tarfile, "-C", exdir))
  if (!identical(status, 0L)) {
    stop("Failed to extract zstd tar archive: ", tarfile)
  }
}


update_matrix_dir <- function(data, old_dir, new_dir) {

  assays <- names(data@assays)
  for (assay.update in assays) {
    data@assays[[assay.update]]@layers <- lapply(
      data@assays[[assay.update]]@layers,
      replace_matrix_dir_paths,
      old_dir = old_dir,
      new_dir = new_dir
    )
  }
  return(data)
}

# Update BPCells matrix paths after extraction
replace_matrix_dir_paths <- function(obj, old_dir, new_dir) {
  new_dir <- normalizePath(new_dir)
  if (!dir.exists(new_dir))
    stop(new_dir, " doesn't exist. Please move BPcells folder first.")

  if (inherits(obj, "MatrixDir") && obj@dir == old_dir) {
    obj@dir <- new_dir
    return(obj)
  }
  if (is.list(obj)) {
    return(lapply(obj, replace_matrix_dir_paths, old_dir, new_dir))
  }
  for (sn in slotNames(obj)) {
    slot(obj, sn) <- replace_matrix_dir_paths(slot(obj, sn), old_dir, new_dir)
  }
  return(obj)
}

# get current matrix_dirs
get_matrix_dirs <- function(scdata) {

  matrix_dirs <- lapply(
    scdata@assays$RNA@layers,
    find_matrix_dir_paths
  )

  matrix_dirs <- unique(unlist(matrix_dirs))
  return(matrix_dirs)
}

find_matrix_dir_paths <- function(obj) {
  matrix_dir_paths <- c()
  if (inherits(obj, "MatrixDir")) {
    return(obj@dir)
  }
  if (is.list(obj)) {
    res <- lapply(obj, find_matrix_dir_paths)
    return(unlist(res))
  }
  for (sn in slotNames(obj)) {
    matrix_dir_paths <- c(
      matrix_dir_paths,
      find_matrix_dir_paths(slot(obj, sn)))
  }
  return(matrix_dir_paths)
}

# Convert disk-backed IterableMatrix to dgCMatrix
#
# Materializes BPCells disk-backed matrices to in-memory dgCMatrix format.
# This is necessary before saving/downloading a Seurat object since the
# disk directories won't be available in the downloaded file.
#
# @param data Seurat object with potential IterableMatrix layers
#
# @return Seurat object with all IterableMatrix converted to dgCMatrix
#
# @export
materialize_bpcells_matrix <- function(data) {
  for (assay_name in Seurat::Assays(data)) {
    assay <- data[[assay_name]]

    # Convert matrix slots
    # Handle layers (counts, data, scale.data, etc.)
    for (layer_name in SeuratObject::Layers(assay)) {
      layer <- SeuratObject::LayerData(assay, layer_name)
      if (is(layer, "IterableMatrix")) {
        SeuratObject::LayerData(data[[assay_name]], layer_name) <-
          as(layer, "dgCMatrix")
      }
    }
  }
  return(data)
}
