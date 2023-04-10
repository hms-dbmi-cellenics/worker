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
  pbmc_small <- SeuratObject::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  pbmc_small@misc$color_pool <- mock_color_pool(20)

  pbmc_small <- Seurat::NormalizeData(pbmc_small, normalization.method = "LogNormalize", verbose = FALSE)
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)

  # scale and PCA
  pbmc_small <- Seurat::NormalizeData(pbmc_small, normalization.method = "LogNormalize", verbose = FALSE)
  pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::RunPCA(pbmc_small, verbose = FALSE, npcs = 10)
  pbmc_small@misc[["active.reduction"]] <- "pca"

  # run UMAP
  pbmc_small <- suppressWarnings(Seurat::RunUMAP(pbmc_small, dims = 1:10, verbose = FALSE))

  return(pbmc_small)
}

mock_cellset_from_python <- function(data) {
  cell_sets <- list(
    "louvain" = list(
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
    ),
    "scratchpad" = list(
      list(
        "key" = "scratchpad-0", "name" = "Custom 0", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#77aadd", "cellIds" = unname(sample(data$cells_id, 5))
      ),
      list(
        "key" = "scratchpad-1", "name" = "Custom 1", "rootNode" = FALSE,
        "type" = "cellSets", "color" = "#ee8866", "cellIds" = unname(sample(data$cells_id, 10))
      )
    ),
    "sample" = list(
      list(
        "key" = "a636ec18-4ba3-475b-989d-0a5b2", "name" = "P13 Acute MISC",
        "color" = "#77aadd", "cellIds" = unname(data$cells_id[1:length(data$cells_id) / 2])
      ),
      list(
        "key" = "eec701f2-5762-4b4f-953d-6aba8", "name" = "P13 Convalescent MISC",
        "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
      )
    ),
    "MISC_status" = list(
      list(
        "key" = "MISC_status-Acute", "name" = "Acute", "color" = "#77aadd",
        "cellIds" = unname(data$cells_id[1:length(data$cells_id) / 2])
      ),
      list(
        "key" = "MISC_status-Convalescent", "name" = "Convalescent",
        "color" = "#ee8866", "cellIds" = unname(data$cells_id[(length(data$cells_id) / 2 + 1):(length(data$cells_id))])
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
    function(i) embedding_data[i,]
  )
  return(embedding_data)
}

mock_req <- function(data) {
  req <- list(
    body = list(
      cellSets = mock_cellset_from_python(data),
      embedding = mock_embedding_data(data)
    )
  )

  return(req)
}

stub_saveRDS <- function(data,fpath) {
  rdir <- file.path(tempdir(), "R")
  if (!dir.exists(rdir)) {
    dir.create(rdir, recursive = TRUE)
  }
  fpath <- file.path(rdir, "r.rds")
  saveRDS(data, fpath)
}

mockery::stub(DownloadAnnotSeuratObject,
              "saveRDS",
              stub_saveRDS)

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
