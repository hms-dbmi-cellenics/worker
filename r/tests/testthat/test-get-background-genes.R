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

  res_gene_counts <- Matrix::rowSums(data@assays$RNA$counts)[res_gene_ids]

  expect_true(all(res_gene_counts > min.total.count))
})

test_that("No genes outside those found have more than min.count counts", {
  data <- mock_scdata()
  req <- mock_req()
  min.total.count <- 15
  res <- getBackgroundExpressedGenes(req, data)
  gene_annotations <- data@misc$gene_annotations

  other_gene_ids <- gene_annotations[!match(res$genes, gene_annotations$name), "input"]

  other_gene_counts <- Matrix::rowSums(data@assays$RNA$counts)[other_gene_ids]

  expect_true(all(other_gene_counts <= min.total.count))
})

test_that("getBackgroundExpressedGenes works with bpcells", {
  data <- mock_scdata(use_bpcells = TRUE)
  req <- mock_req()
  min.total.count <- 15
  res <- getBackgroundExpressedGenes(req, data)
  gene_annotations <- data@misc$gene_annotations

  other_gene_ids <- gene_annotations[!match(res$genes, gene_annotations$name), "input"]

  other_gene_counts <- Matrix::rowSums(data@assays$RNA$counts)[other_gene_ids]

  expect_true(all(other_gene_counts <= min.total.count))
})