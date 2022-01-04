mock_req <- function(type = "louvain") {
  req <- list(
    body =
      list(
        config = list(
          resolution = 2
        ),
        type = type
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
  return(pbmc_small)
}

test_that("clustering returns correct keys", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_keys <- c("cluster", "cell_ids")

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- runClusters(req, data)
    expect_equal(names(res), expected_keys)
  }
})

test_that("clustering returns one value per cell", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_n_cells <- ncol(data)

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- runClusters(req, data)
    n_cells <- nrow(res)

    expect_equal(n_cells, expected_n_cells)
  }
})

test_that("clustering orders barcodes correctly", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_barcodes <- colnames(data)

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- runClusters(req, data)
    barcodes <- rownames(res)
    expect_equal(barcodes, expected_barcodes)
  }
})

test_that("clustering returns at least one cluster", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- runClusters(req, data)
    n_clusters <- length(unique(res$cluster))
    expect_gte(n_clusters, 1)
  }
})
