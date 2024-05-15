mock_scdata <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  gene_annotations <- data.frame(
    input = paste0("ENSG", seq_len(nrow(pbmc_raw))),
    name = row.names(pbmc_raw),
    row.names = paste0("ENSG", seq_len(nrow(pbmc_raw)))
  )
  gene_annotations$original_name <- gene_annotations$name
  row.names(pbmc_raw) <- gene_annotations$input

  pbmc_raw <- as(as.matrix(pbmc_raw), 'dgCMatrix')
  pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  return(pbmc_small)
}

mock_req <- function() {
  backgroundCells <- 0:39
  baseCells <- 40:79

  req <- list(body = list(backgroundCells = backgroundCells, baseCells = baseCells))
  return(req)
}

test_that("addComparisonGroup properly identifies the cells", {
  data <- mock_scdata()
  req <- mock_req()

  data <- addComparisonGroup(req, data)
  expect_true(all(data$custom[1:40] == "background"))
  expect_true(all(data$custom[41:80] == "base"))
})

test_that("All genes found in getBackgroundExpressedGenes have more than min.count counts", {
  data <- mock_scdata()
  req <- mock_req()
  min.total.count <- 15
  res <- getBackgroundExpressedGenes(req, data)
  gene_annotations <- data@misc$gene_annotations

  res_gene_ids <- gene_annotations[match(res$genes, gene_annotations$name), "input"]

  res_gene_counts <- Matrix::rowSums(data@assays$RNA@counts)[res_gene_ids]

  expect_true(all(res_gene_counts > min.total.count))
})

test_that("No genes outside those found have more than min.count counts", {
  data <- mock_scdata()
  req <- mock_req()
  min.total.count <- 15
  res <- getBackgroundExpressedGenes(req, data)
  gene_annotations <- data@misc$gene_annotations

  other_gene_ids <- gene_annotations[!match(res$genes, gene_annotations$name), "input"]

  other_gene_counts <- Matrix::rowSums(data@assays$RNA@counts)[other_gene_ids]

  expect_true(all(other_gene_counts <= min.total.count))
})
