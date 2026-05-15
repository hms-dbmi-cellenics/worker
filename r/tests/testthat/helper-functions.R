mock_color_pool <- function(n) {
  paste0("color_", 1:n)
}

mock_scdata <- function(
  counts = NULL,
  with_umap = FALSE,
  use_bpcells = FALSE,
  use_enids = TRUE,
  nreps = 1
) {

  if (is.null(counts)) {
    counts <- read.table(
      file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
      as.is = TRUE
    )
  }

  if (nreps > 1) {
    counts <- do.call("cbind", replicate(nreps, counts, simplify = FALSE))
    colnames(counts) <- make.unique(colnames(counts))
  }

  # some tests use default rownames for simplicity
  if (use_enids) {
    input_names <- paste0("ENSG", seq_len(nrow(counts)))
  } else {
    input_names <- row.names(counts)
  }

  gene_annotations <- data.frame(
    input = input_names,
    name = row.names(counts),
    original_name = row.names(counts),
    row.names = input_names
  )

  row.names(counts) <- input_names
  counts <- as(as.matrix(counts), "dgCMatrix")

  if (use_bpcells) {
    matrix_dir <- file.path(tempdir(), "matrix_dir")
    if (dir.exists(matrix_dir)) {
      unlink(matrix_dir, recursive = TRUE)
    }
    counts <- BPCells::convert_matrix_type(counts)
    counts <- BPCells::write_matrix_dir(counts, matrix_dir)
  }

  scdata <- SeuratObject::CreateSeuratObject(counts = counts)

  scdata$cells_id <- 0:(ncol(scdata) - 1)
  scdata@misc$gene_annotations <- gene_annotations
  scdata@misc$color_pool <- mock_color_pool(20)

  scdata <- Seurat::NormalizeData(
    scdata,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
  scdata <- Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)

  scdata <- Seurat::RunPCA(scdata, verbose = FALSE, npcs = 10)
  scdata@misc[["active.reduction"]] <- "pca"
  scdata@misc[["numPCs"]] <- 10

  if (with_umap) {
    scdata <- suppressWarnings(
      Seurat::RunUMAP(scdata, dims = 1:10, verbose = FALSE)
    )
  }

  # add gene dispersions
  vars <- Seurat::HVFInfo(object = scdata, assay = "RNA", method = "vst")
  vars$SYMBOL <- gene_annotations$name[
    match(rownames(vars),
    gene_annotations$input)
  ]
  vars$ENSEMBL <- rownames(vars)

  scdata@misc[["gene_dispersion"]] <- vars

  # add mitochondrial content and doublet scores
  set.seed(0)
  scdata[["percent.mt"]] <- rnorm(ncol(scdata), 5, 1)
  scdata[["doublet_scores"]] <- rnorm(ncol(scdata), 0.5, 0.1)

  return(scdata)
}

mock_sketch <- function(scdata, ncells = 100) {
  set.seed(123)
  sampled_cells <- sample(Seurat::Cells(scdata), size = ncells)
  sketch_assay <- subset(
    scdata[[Seurat::DefaultAssay(scdata)]],
    cells = sampled_cells
  )
  scdata[["sketch"]] <- sketch_assay
  Seurat::DefaultAssay(scdata) <- "sketch"
  return(scdata)
}
