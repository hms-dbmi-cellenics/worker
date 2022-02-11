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

  row.names(pbmc_raw) <- gene_annotations$input


  pbmc_small <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  pbmc_small@misc$color_pool <- list("#e377c2","#8c564b","#d62728","#2ca02c","#ff7f0e")
  return(pbmc_small)
}

test_that("getExpressionCellSetIDs with a single expression filter returns correct cell ids", {
  data <- mock_scdata()
  req <- list(list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5))
  res <- getExpressionCellSetIDs(req, data)

  annot <- data@misc$gene_annotations
  enid <- annot$input[annot$name == req[[1]]$geneName]
  expected <- data$cells_id[data[["RNA"]]@data[enid, ] > 0.5]
  expect_equal(res$keep_ids, expected)
})

test_that("getExpressionCellSetIDs with two expression filters returns correct cell ids", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0.5)
  )
  res <- getExpressionCellSetIDs(req, data)

  # gives same results as manual
  annot <- data@misc$gene_annotations
  enid1 <- annot$input[annot$name == req[[1]]$geneName]
  enid2 <- annot$input[annot$name == req[[2]]$geneName]
  expected <- data$cells_id[data[["RNA"]]@data[enid1, ] > 0.5 & data[["RNA"]]@data[enid2, ] > 0.5]

  expect_equal(res$keep_ids, expected)

  # order invariant
  req2 <- list(
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5)
  )

  res2 <- getExpressionCellSetIDs(req2, data)
  expect_equal(res$keep_ids, res2$keep_ids)
})


test_that("getExpressionCellSetIDs returns no cell ids if filters remove everything", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "MS4A1", comparisonType = "lessThan", thresholdValue = 0.5)
  )
  res <- getExpressionCellSetIDs(req, data)
  expect_length(res$keep_ids, 0)
})


test_that("getExpressionCellSetIDs returns all cell ids if filters include everything", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = -1)
  )
  res <- getExpressionCellSetIDs(req, data)
  expect_equal(res$keep_ids, data$cells_id)
})



test_that("getExpressionCellSetIDs fails if requested geneNames are not present in SeuratObject", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "BLAH", comparisonType = "greaterThan", thresholdValue = 0)
  )

  expect_error(getExpressionCellSetIDs(req, data), "gene name\\(s\\) that are not present[.]")
})

test_that("CellSet naming is correct", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0),
    list(geneName = "MS4A1", comparisonType = "lessThan", thresholdValue = 0.5)
  )

  res <- getExpressionCellSetIDs(req, data)

  expect_equal(res$cell_set_name, "CD79B>0, MS4A1<0.5")
})



