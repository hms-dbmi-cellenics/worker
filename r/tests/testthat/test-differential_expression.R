mock_req <- function() {
  req <- list(body = list(
    comparisonType = "within",
    backgroundCells = 0:39,
    baseCells = 40:79,
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
    comparisonType = "within",
    backgroundCells = 0:39,
    baseCells = 40:79,
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
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  enids <- paste0("ENSG", seq_len(nrow(pbmc_raw)))
  gene_annotations <- data.frame(
    input = enids,
    name = row.names(pbmc_raw),
    row.names = enids
  )

  row.names(pbmc_raw) <- enids
  pbmc_small <- SeuratObject::CreateSeuratObject(counts = pbmc_raw)

  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  return(pbmc_small)
}

test_that("runDE generates the expected return format for comparisons within samples/groups", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDE(req, data)
  res$gene_results <- purrr::transpose(res$gene_results)

  # number of genes is number of possible DE rows
  expect_equal(res$full_count, nrow(data))

  # returning only at most limit number of genes
  expect_equal(length(res$gene_results), req$body$pagination$limit)

  # ordering is correct
  logfc <- sapply(res$gene_results, `[[`, "logFC")
  expect_equal(logfc, sort(logfc, decreasing = TRUE))

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

test_that("runDE generates the expected return format for comparisons between samples/groups with 3+ samples", {
  data <- mock_scdata()
  data$samples <- rep(LETTERS[1:4], each = 20)
  req <- mock_req()
  req$body$comparisonType <- "between"



  res <- runDE(req, data)
  res$gene_results <- purrr::transpose(res$gene_results)

  # number of genes is the min of pagination limit and DE genes
  expect_equal(res$full_count, min(nrow(data), req$body$pagination$limit))

  # returning only at most limit number of genes
  expect_equal(length(res$gene_results), req$body$pagination$limit)

  # ordering is correct
  logfc <- sapply(res$gene_results, `[[`, "logFC")
  expect_equal(logfc, sort(logfc, decreasing = TRUE))

  # have the correct column names
  expect_columns <-
    c(
      "p_val",
      "logFC",
      "AveExpr",
      "p_val_adj",
      "gene_names",
      "Gene"
    )
  expect_equal(unique(names(unlist(res$gene_results))), expect_columns)
})

test_that("runDE generates the expected return format for comparisons between samples/groups with less than 3 samples", {
  data <- mock_scdata()
  data$samples <- rep(LETTERS[1:2], each = 40)
  req <- mock_req()
  req$body$comparisonType <- "between"



  res <- runDE(req, data)
  res$gene_results <- purrr::transpose(res$gene_results)

  # number of genes is the min of pagination limit and DE genes
  expect_equal(res$full_count, min(nrow(data), req$body$pagination$limit))

  # returning only at most limit number of genes
  expect_equal(length(res$gene_results), req$body$pagination$limit)

  # ordering is correct
  logfc <- sapply(res$gene_results, `[[`, "logFC")
  expect_equal(logfc, sort(logfc, decreasing = TRUE))

  # have the correct column names
  expect_columns <-
    c(
      "logFC",
      "AveExpr",
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
  res$gene_results <- purrr::transpose(res$gene_results)
  expect_equal(length(res$gene_results), nrow(data))
})

test_that("runDE was able to convert from ensemblIDs to gene symbols", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDE(req, data)
  res$gene_results <- purrr::transpose(res$gene_results)

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
  expect_equal(res$gene_results$gene_names[[1]], "CST3")
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
  expect_equal(length(res$gene_results$Gene), 0)

  req$body$pagination$filters <-
    list(list(
      columnName = "logFC",
      comparison = "greaterThan",
      value = Inf
    ))

  res <- runDE(req, data)
  expect_equal(length(res$gene_results$Gene), 0)
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
