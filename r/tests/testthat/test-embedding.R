library('mockery')

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

test_that("TSNE embedding works", {

  mock_RunTSNE <- function(config, method, reduction_type, num_pcs, data) {
    Seurat::RunTSNE(data,
     reduction = reduction_type,
     dims = 1:num_pcs,
     perplexity = config$perplexity,
     learning.rate = config$learningRate,
     check_duplicates = FALSE
    )
  }

  reduction_method = 'tsne'
  stub(runEmbedding, 'getEmbedding', mock_RunTSNE)

  data <- suppressWarnings(mock_scdata())
  req <- list(
    body = list(
      type = reduction_method,
      config = list(perplexity = 10, learningRate = 100)
    )
  )

  res <- runEmbedding(req, data)
  expected_res <- as.data.frame(Seurat::Embeddings(data)[,1:2])

  # Expect all cells to be in embedding
  expect_equal(length(res), length(expected_res$PC_1))
})

test_that("UMAP embedding works", {
  mock_RunUMAP <- function(config, method, reduction_type, num_pcs, data) {
    Seurat::RunUMAP(data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      verbose = FALSE,
      min.dist = config$minimumDistance,
      metric = config$distanceMetric,
      seed.use = ULTIMATE_SEED
    )
  }

  stub(runEmbedding, 'getEmbedding', mock_RunUMAP)

  reduction_method <- "umap"

  data <- suppressWarnings(mock_scdata())
  req <- list(
    body = list(
      type = reduction_method,
      config = list(minimumDistance = 0.1, distanceMetric = "cosine")
    )
  )

  res <- runEmbedding(req, data)

  expected_res <- as.data.frame(Seurat::Embeddings(data)[,1:2])

  # Expect all cells to be in embedding
  expect_equal(length(res), length(expected_res$PC_1))
})

test_that("RunTSNE uses the correct params", {

  mock_RunTSNE <- mock(TRUE)

  stub(getEmbedding, 'Seurat::RunTSNE', mock_RunTSNE)

  data <- suppressWarnings(mock_scdata())
  config <- list(perplexity = 10, learningRate = 100)
  reduction_type <- "pca"
  method <- "tsne"
  num_pcs <- 1

  res <- getEmbedding(config, method, reduction_type, num_pcs, data)

  # Check that tsne is called using the correct parameters
  expect_equal(length(mock_RunTSNE), 1)
  args <- mock_args(mock_RunTSNE)

  expect_equal(args[[1]], list(data,
                               reduction = reduction_type,
                               dims = 1:num_pcs,
                               perplexity = config$perplexity,
                               learning.rate = config$learningRate))

})

test_that("RunUMAP uses umap-learn with seed.use", {

  mock_RunUMAP <- mock(TRUE)

  stub(getEmbedding, 'Seurat::RunUMAP', mock_RunUMAP)

  data <- suppressWarnings(mock_scdata())
  config <- list(minimumDistance = 0.1, distanceMetric = "cosine")
  reduction_type <- "pca"
  method <- "umap"
  num_pcs <- 1

  res <- getEmbedding(config, method, reduction_type, num_pcs, data)

  # Check that umap is called using the correct parameters
  expect_equal(length(mock_RunUMAP), 1)
  args <- mock_args(mock_RunUMAP)

  expect_equal(args[[1]], list(data,
   reduction = reduction_type,
   dims = 1:num_pcs,
   verbose = FALSE,
   min.dist = config$minimumDistance,
   metric = config$distanceMetric,
   umap.method = "umap-learn",
   seed.use = ULTIMATE_SEED))

})
