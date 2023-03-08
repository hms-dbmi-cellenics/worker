mock_gene_results <- function() {
  data(exprs, package = "presto")
  data(y, package = "presto")

  result <- presto::wilcoxauc(exprs, y)
  result <- result[result$group == "A", ]
  rownames(result) <- result$feature
  result$gene_names <- result$feature
  result <- result[, c("gene_names", "pval", "logFC", "pct_in", "pct_out", "padj", "auc")]
  colnames(result) <- list(
    "gene_names", "p_val", "logFC", "pct_1", "pct_2", "p_val_adj", "auc"
  )
  return(result)
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )

  # change a few names, to ease testing getMarkerNames

  ngenes <- nrow(pbmc_small)
  sampled_genes <- seq(1, ngenes, by = 5)
  new_names <- paste0("name_", sampled_genes)
  pbmc_small@misc$gene_annotations$name[sampled_genes] <- new_names

  return(pbmc_small)
}

mock_cellSets <- function() {
  cellSets <- list(
    children = list(
      louvain1 = list(cellIds = 0:39),
      louvain2 = list(cellIds = 40:79)
    )
  )
  cellSets
}

test_that("applyFilters works with gene with exact match", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "G5"))

  gene_results <- applyFilters(gene_results, filters)

  expect_equal(gene_results$gene_names, "G5")
})

test_that("applyFilters works with gene with case insensitive match", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "g5"))

  gene_results <- applyFilters(gene_results, filters)

  expect_equal(gene_results$gene_names, "G5")
})

test_that("applyFilters returns all results that match gene search pattern", {
  gene_results <- mock_gene_results()
  filters <- list(list(columnName = "gene_names", expression = "G"))

  filt_results <- applyFilters(gene_results, filters)

  expect_equal(filt_results, gene_results)
})


test_that("applyFilters works with single numeric filter", {
  gene_results <- mock_gene_results()
  filter_up <- list(list(
    columnName = "logFC",
    comparison = "greaterThan", value = 0
  ))
  filter_dn <- list(list(
    columnName = "logFC",
    comparison = "lessThan", value = 0
  ))

  results_up <- applyFilters(gene_results, filter_up)
  results_dn <- applyFilters(gene_results, filter_dn)

  expect_equal(results_up, gene_results[gene_results$logFC > 0, ])
  expect_equal(results_dn, gene_results[gene_results$logFC < 0, ])
})

test_that("running applyFilters sequentially is equivalent to running together", {
  gene_results <- mock_gene_results()
  filter_gene <- list(columnName = "gene_names", expression = "G1")
  filter_up <- list(columnName = "logFC", comparison = "greaterThan", value = 0)

  # apply sequentially
  results_seq <- applyFilters(gene_results, list(filter_gene))
  results_seq <- applyFilters(results_seq, list(filter_up))

  # apply together
  results_together <- applyFilters(gene_results, list(filter_gene, filter_up))

  expect_equal(results_seq, results_together)
})

# subsetIds

test_that("subsetIds returns same class as input data", {
  data <- mock_scdata()
  cell_ids <- c(1, 2, 3, 5, 7, 11, 13, 17, 19)

  data_subset <- subsetIds(data, cell_ids)

  expect_true(class(data) == class(data_subset))
})

test_that("subsetIds returns seurat object correctly", {
  data <- mock_scdata()
  cell_ids <- c(1, 2, 3, 5, 7, 11, 13, 17, 19)

  data_subset <- subsetIds(data, cell_ids)

  expect_equal(length(cell_ids), ncol(data_subset))
  expect_equal(data_subset@meta.data$cells_id, cell_ids)
  expect_true(all(colnames(data_subset) %in% colnames(data)))
})

test_that("subsetIds errors when empty cell_ids", {
  data <- mock_scdata()
  cell_ids <- c()

  expect_error(subsetIds(data, cell_ids))
})

# getTopMarkerGenes

test_that("getTopMarkerGenes returns an object with correct columns", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 5

  expected_columns <-
    c(
      "feature",
      "group",
      "avgExpr",
      "logFC",
      "statistic",
      "auc",
      "pval",
      "padj",
      "pct_in",
      "pct_out"
    )

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  res <- getTopMarkerGenes(nFeatures, data, cell_sets_ids)

  expect_equal(colnames(res), expected_columns)
})

test_that("getTopMarkerGenes returns at least 1 gene and n_genes * n_cellSets at the most", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 1000

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  res <- getTopMarkerGenes(nFeatures, data, cell_sets_ids)

  expect_gte(nrow(res), 1)
  expect_lte(nrow(res), nFeatures * length(cellSets))
})

test_that("getTopMarkerGenes returns correctly filtered marker genes", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  aucMin <- 0.2
  pctInMin <- 17
  pctOutMax <- 73

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  res <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids,
      aucMin,
      pctInMin,
      pctOutMax
    )

  expect_true(all(res$auc >= aucMin))
  expect_true(all(res$pct_in >= pctInMin))
  expect_true(all(res$pct_out <= pctOutMax))
})

test_that("getTopMarkerGenes returns empty if no genes match filters", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  aucMin <- 1
  pctInMin <- 100
  pctOutMax <- 0

  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  res <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids,
      aucMin,
      pctInMin,
      pctOutMax
    )

  expect_equal(nrow(res), 0)
})

test_that("getTopMarkerGenes markers in correct order", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  
  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  res <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids
    )

  expect_equal(res, dplyr::arrange(res, group))
})

# getMarkerNames

test_that("getMarkerNames returns object with correct columns", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  all_markers <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids
    )

  res <- getMarkerNames(data, all_markers)
  expected_columns <- c("group", "input", "name")

  expect_equal(names(res), expected_columns)
})

test_that("getMarkerNames returns same number of gene names as requested", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  all_markers <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids
    )

  res <- getMarkerNames(data, all_markers)

  expect_equal(nrow(res), nrow(all_markers))
})

test_that("getMarkerNames correct gene names", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  all_markers <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids
    )

  res <- getMarkerNames(data, all_markers)

  expected_names <- data@misc$gene_annotations %>%
    dplyr::filter(input %in% all_markers$feature) %>%
    dplyr::pull(name)

  # getTopMarkerGenes orders them by group
  expect_setequal(res$name, expected_names)
})

test_that("getMarkerNames returns input if there's no gene name", {
  data <- mock_scdata()
  cellSets <- mock_cellSets()$children
  nFeatures <- 42
  cell_sets_ids <- lapply(cellSets, function(x) x[["cellIds"]])

  all_markers <-
    getTopMarkerGenes(
      nFeatures,
      data,
      cell_sets_ids
    )

  res <- getMarkerNames(data, all_markers)

  expected_noname_genes <- data@misc$gene_annotations %>%
    dplyr::filter(!grepl("^name_", name))

  named_genes <- res %>%
    dplyr::filter(grepl("^name_", name))

  noname_genes <- res %>%
    dplyr::filter(!grepl("^name_", name))

  expect_true(all(noname_genes$input %in% expected_noname_genes$input))
  expect_false(any(named_genes$input %in% expected_noname_genes$input))
})


# getSNNigraph

# not sure how to test this
test_that("getSNNigraph returns an igraph object with correct dimensions", {
  data <- mock_scdata()
  g <- getSNNiGraph(data)
})

test_that("generateErrorMessage concatenates error code and user message with :|:", {
  code <- "error_code"
  msg <- "user_message"
  expect_equal(generateErrorMessage(code, msg), "error_code:|:user_message")
})

test_that("extractErrorList extracts error code and user message seperated by :|:", {
  error_list <- extractErrorList("error_code:|:user_message")
  expect_equal(error_list$user_message, "user_message")
  expect_equal(error_list$error_code, "error_code")
})

test_that("extractErrorList correctly formats unhandled error messages", {
  unhandled_message <- tryCatch(stop("unhandled message!"), error = function(e) {
    return(e$message)
  })

  error_list <- extractErrorList(unhandled_message)
  expect_equal(error_list$user_message, unhandled_message)
  expect_equal(error_list$error_code, "R_WORKER_UNHANDLED_ERROR")
})

test_that("formatResponse creates a list with data and error objects", {
  expect_equal(
    formatResponse(letters[1:5], "error!"),
    list(
      data = letters[1:5],
      error = "error!"
    )
  )
})
