#' Annotate Seurat object
#'
#' This function adds information about embeddings and clustering to the processed Seurat object.
#' Embedding are obtained and added to the Seurat object using the assignEmbedding function.
#' Cellsets object is received from the Python worker and parsed to extract
#' relevant information about louvain/leiden clusters, custom cellsets,
#' and cluster annotations. This information is added to the Seurat object
#' to provide a more comprehensive representation of the data.
#'
#'
#' @param req Request object that contains the cellsets object and embedding data.
#' @param data Seurat object
#'
#' @return Path to the saved zip with the Seurat object.
#' @export
#'
DownloadAnnotSeuratObject <- function(req, data) {
  cell_sets <- req$body$cellSets

  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  if (!req$body$isObj2s) {
    embedding_data <- req$body$embedding
    embedding_method <- req$body$embeddingMethod
    data <- assignEmbedding(embedding_data, data, embedding_method)
  }

  # disk-backed (BPCells) counts stay on disk: reading them into memory crashes
  # the worker for large experiments, so their matrix dir is shipped alongside
  save_download_zip(data, DOWNLOAD_ZIP_PATH)
  return(DOWNLOAD_ZIP_PATH)
}


#' Save a Seurat object for download as a zip
#'
#' The zip holds \code{object.rds}, the BPCells matrix dirs its disk-backed
#' layers read from (\code{matrix_dir}) and \code{load_object.R}, which loads the
#' object and points those layers at the extracted matrix dirs.
#'
#' @param data Seurat object
#' @param zip_path path of the zip to write
#'
#' @return \code{zip_path}
#' @export
#'
save_download_zip <- function(data, zip_path) {
  out_dir <- tempfile("download_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # copy each matrix dir next to the object, as matrix_dir(_2, ...)
  matrix_dirs <- get_matrix_dirs(data)
  for (i in seq_along(matrix_dirs)) {
    name <- if (i == 1) "matrix_dir" else paste0("matrix_dir_", i)
    new_dir <- file.path(out_dir, name)
    dir.create(new_dir)
    file.copy(list.files(matrix_dirs[i], full.names = TRUE), new_dir, recursive = TRUE)
    data <- update_matrix_dir(data, matrix_dirs[i], new_dir)
  }

  saveRDS(data, file.path(out_dir, "object.rds"))
  writeLines(DOWNLOAD_LOAD_SCRIPT, file.path(out_dir, "load_object.R"))

  # stored, not compressed: matrix dirs and the rds are compressed already
  if (file.exists(zip_path)) unlink(zip_path)
  zip_path <- normalizePath(zip_path, mustWork = FALSE)
  workdir <- setwd(out_dir)
  on.exit(setwd(workdir), add = TRUE)
  utils::zip(zip_path, files = list.files("."), flags = "-r0q")

  return(zip_path)
}


DOWNLOAD_LOAD_SCRIPT <- '# Load the Seurat object downloaded from Cellenics.
#
# Its counts are disk-backed (BPCells): they stay in matrix_dir and are read
# from there as needed, so large experiments fit in memory. Keep matrix_dir next
# to object.rds, and install BPCells (https://bnprks.github.io/BPCells):
#   remotes::install_github("bnprks/BPCells/r")

library(Seurat)
library(BPCells)

# folder this zip was extracted to
data_dir <- "."

scdata <- readRDS(file.path(data_dir, "object.rds"))

# point the disk-backed layers at the matrix dirs in data_dir
update_matrix_dirs <- function(obj) {
  if (inherits(obj, "MatrixDir")) {
    obj@dir <- normalizePath(file.path(data_dir, basename(obj@dir)))
    return(obj)
  }
  if (is.list(obj)) return(lapply(obj, update_matrix_dirs))
  for (slot_name in methods::slotNames(obj)) {
    methods::slot(obj, slot_name) <- update_matrix_dirs(methods::slot(obj, slot_name))
  }
  obj
}

for (assay in names(scdata@assays)) {
  scdata@assays[[assay]]@layers <- lapply(scdata@assays[[assay]]@layers, update_matrix_dirs)
}

# To read the counts into memory instead (needs enough RAM for them), e.g. to
# re-upload the object to Cellenics:
# scdata[["RNA"]]$counts <- as(scdata[["RNA"]]$counts, "dgCMatrix")
'
