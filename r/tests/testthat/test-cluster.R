mock_req <- function(type = "louvain") {
  req <- list(
    body =
      list(
        config = list(
          resolution = 0.5
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

test_that("leiden clustering works", {

})
