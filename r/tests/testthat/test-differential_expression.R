mock_req <- function() {
  req <- list(body = list(
    backgroundCells = 0:50,
    baseCells = 51:100,
    pagination = list(
      orderBy = "logFC",
      orderDirection = "DESC",
      offset = 0,
      limit = 50
    )
  ))
}

mock_req_genes_only <- function() {
  req <- list(body = list(
    backgroundCells = 0:50,
    baseCells = 51:100,
    pagination = list(
      orderBy = "logFC",
      orderDirection = "DESC",
      offset = 0,
      limit = 200
    ),
    genesOnly = TRUE
  ))
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

test_that("runDE generates the expected return format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDE(req, data)

  # number of genes is number of possible DE rows
  expect_equal(res$full_count, nrow(data))

  # returning only at most limit number of genes
  expect_equal(length(res$gene_results), req$body$pagination$limit)

  # ordering is correct
  expect_equal(
    res$gene_results$logFC,
    sort(res$gene_results$logFC, decreasing = TRUE)
  )

  # have the correct column names
  expect_columns <-
    c(
      "p_val",
      "logFC",
      "pct_1",
      "pct_2",
      "p_val_adj",
      "auc",
      "gene_names",
      "Gene"
    )
  expect_equal(unique(names(unlist(res$gene_results))), expect_columns)
})

test_that("runDE limit won't return more than available genes", {
  data <- mock_scdata()
  req <- mock_req()

  req$body$pagination$limit <- nrow(data) + 50

  res <- runDE(req, data)
  expect_equal(length(res$gene_results), nrow(data))
})

test_that("runDE was able to convert from ensemblIDs to gene symbols", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDE(req, data)
  expect_equal(sum(is.na(res$gene_results$gene_names)), 0)

  # runDE returns list of named lists
  res_df <- dplyr::bind_rows(res$gene_results)

  # check it has at least one gene (previously got NULL and next test passed)
  expect_type(res_df$gene_names, "character")
  expect_gte(length(res_df$gene_names), 1)

  expect_true(all(res_df$gene_names %in% data@misc$gene_annotations$name))
  expect_equal(length(unique(res_df$gene_names)), length(res$gene_results))
})

test_that("runDE works with gene name filter", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$pagination$filters <-
    list(list(columnName = "gene_names", expression = "CST3"))


  res <- runDE(req, data)
  expect_equal(res$gene_results[[1]]$gene_names, "CST3")
})

test_that("runDE works with numeric filters", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$pagination$filters <-
    list(list(
      columnName = "logFC",
      comparison = "greaterThan",
      value = 0
    ))


  res <- runDE(req, data)
  res_df <- dplyr::bind_rows(res$gene_results)
  expect_true(all(res_df$logFC > 0))
})

test_that("runDE works when no genes match filters", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$pagination$filters <-
    list(list(columnName = "gene_names", expression = "mythical_gene"))

  res <- runDE(req, data)
  expect_equal(length(res$gene_results), 0)

  req$body$pagination$filters <-
    list(list(
      columnName = "logFC",
      comparison = "greaterThan",
      value = Inf
    ))

  res <- runDE(req, data)
  expect_equal(length(res$gene_results), 0)
})

test_that("DE with genes_only returns list of ENSEMBLIDS and gene symbols ", {
  data <- mock_scdata()
  req <- mock_req_genes_only()

  res <- runDE(req, data)

  expect_equal(length(res), 2)
  expect_true("gene_results" %in% names(res))

  gene_results_names <- names(res$gene_results)
  expect_true("gene_id" %in% gene_results_names & "gene_names" %in% gene_results_names)

  expect_equal(length(res$gene_results[[1]]), res$full_count)
})
