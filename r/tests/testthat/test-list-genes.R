mock_req <- function(orderBy = "dispersions", orderDirection = "DESC", offset = 0, limit = 20, geneNamesFilter = NULL) {
  req <- list(
    body = list(
      "name" = "ListGenes",
      "selectFields" = list("gene_names", "dispersions"),
      "orderBy" = orderBy,
      "orderDirection" = orderDirection,
      "offset" = offset,
      "limit" = limit,
      "geneNamesFilter" = geneNamesFilter
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

  vars <- Seurat::HVFInfo(object = pbmc_small, assay = "RNA", selection.method = "vst")
  annotations <- pbmc_small@misc[["gene_annotations"]]
  vars$SYMBOL <- rownames(vars)
  vars$ENSEMBL <- annotations$input[match(rownames(vars), annotations$name)]
  pbmc_small@misc[["gene_dispersion"]] <- vars

  return(pbmc_small)
}

test_that("List of genes generates the expected return format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getList(req, data)

  # number of genes is number of possible DE rows
  expect_equal(res$full_count, 50)

  # returning only at most limit number of genes
  expect_equal(nrow(res$gene_results), req$body$limit)

  # ordering is correct
  expect_equal(res$gene_results$dispersions, sort(res$gene_results$dispersions, decreasing = TRUE))

  # have the correct column names
  expect_columns <- c("gene_names", "dispersions")
  expect_equal(colnames(res$gene_results), expect_columns)
})

test_that("Order direction works properly", {
  data <- mock_scdata()
  req <- mock_req(orderBy = "dispersions", orderDirection = "DESC", offset = 0, limit = 20)

  res <- getList(req, data)

  # ordering is correct
  expect_equal(res$gene_results$dispersions, sort(res$gene_results$dispersions, decreasing = TRUE))

  req <- mock_req(orderBy = "dispersions", orderDirection = "ASC", offset = 0, limit = 20)

  res <- getList(req, data)

  # ordering is correct
  expect_equal(res$gene_results$dispersions, sort(res$gene_results$dispersions, decreasing = FALSE))
})

test_that("Order by works properly", {
  data <- mock_scdata()
  req <- mock_req(orderBy = "gene_names", orderDirection = "DESC", offset = 0, limit = 20)

  res <- getList(req, data)

  # ordering is correct
  expect_equal(res$gene_results$gene_names, sort(res$gene_results$gene_names, decreasing = TRUE))

  req <- mock_req(orderBy = "gene_names", orderDirection = "ASC", offset = 0, limit = 20)

  res <- getList(req, data)

  # ordering is correct
  expect_equal(res$gene_results$gene_names, sort(res$gene_results$gene_names, decreasing = FALSE))
})

test_that("Pagination works", {
  data <- mock_scdata()
  req <- mock_req(orderBy = "gene_names", orderDirection = "DESC", offset = 20, limit = 40)

  res <- getList(req, data)

  expect_equal(length(res$gene_results$gene_names), 40)

  expect_equal(res$gene_results$gene_names, sort(rownames(data), decreasing = TRUE)[21:60])

  req <- mock_req(orderBy = "dispersions", orderDirection = "DESC", offset = 20, limit = 40)

  res <- getList(req, data)

  expect_equal(length(res$gene_results$gene_names), 40)

  expect_equal(res$gene_results$dispersions, sort(data@misc$gene_dispersion$variance.standardized, decreasing = TRUE)[21:60])
})

test_that("Start with pattern is applied", {
  pat <- "^GZ"
  data <- mock_scdata()
  req <- mock_req(
    orderBy = "gene_names",
    orderDirection = "DESC",
    offset = 0,
    limit = 40,
    geneNamesFilter = pat
  )

  res <- getList(req, data)

  expect_true(all(grepl(pat, res$gene_results$gene_names)))

  grep_results <- grepl(pat, rownames(data))
  expect_true(all(res$gene_results$gene_names %in% data@misc$gene_annotations[grep_results, "name"]))
  expect_equal(res$full_count, sum(grep_results == TRUE))
})

test_that("Ends with pattern is applied", {
  pat <- "1$"
  data <- mock_scdata()
  req <- mock_req(
    orderBy = "gene_names",
    orderDirection = "DESC",
    offset = 0,
    limit = 40,
    geneNamesFilter = pat
  )

  res <- getList(req, data)

  expect_true(all(grepl(pat, res$gene_results$gene_names)))

  grep_results <- grepl(pat, rownames(data))
  expect_true(all(res$gene_results$gene_names %in% data@misc$gene_annotations[grep_results, "name"]))
  expect_equal(res$full_count, 20)
})

test_that("Contains pattern is applied", {
  pat <- "CR"
  data <- mock_scdata()
  req <- mock_req(
    orderBy = "gene_names",
    orderDirection = "DESC",
    offset = 0,
    limit = 40,
    geneNamesFilter = pat
  )

  res <- getList(req, data)

  expect_true(all(grepl(pat, res$gene_results$gene_names)))

  grep_results <- grepl(pat, rownames(data))
  expect_true(all(res$gene_results$gene_names %in% data@misc$gene_annotations[grep_results, "name"]))
  expect_equal(res$full_count, 40)
})
