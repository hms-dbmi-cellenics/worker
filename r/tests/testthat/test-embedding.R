mock_req <- function() {
  req <- list(
    body = list(
      type = "umap",
      config = list(minimumDistance = 0.1, distanceMetric = "cosine"),
      use_saved = FALSE
    )
  )
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = paste0("ENSG", seq_len(nrow(pbmc_small))),
    name = row.names(pbmc_small),
    row.names = paste0("ENSG", seq_len(nrow(pbmc_small)))
  )

  pbmc_small <- Seurat::RunPCA(pbmc_small, npcs = 5)
  pbmc_small$samples <- rep("Sample1", 80)
  pbmc_small@misc$numPCs <- 5
  return(pbmc_small)
}

test_that("PCA embedding works", {
  data <- suppressWarnings(mock_scdata())
  req <- mock_req()
  req$body$type <- "pca"

  res <- runEmbedding(req, data)

  expected_res <- as.data.frame(Seurat::Embeddings(data)[,1:2])

  expected_res <- expected_res %>%
    as.data.frame() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(PCS = list(c(PC_1, PC_2)))

  expect_equal(res,expected_res$PCS)
})
