mock_req <- function() {
  req <- list()
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  # create randomized cell ids (as done in the pipeline)
  set.seed(1)
  pbmc_small$cells_id <- sample(0:(ncol(pbmc_small) - 1))
  pbmc_small@misc$gene_annotations <- data.frame(
    input = paste0("ENSG", seq_len(nrow(pbmc_small))),
    name = row.names(pbmc_small),
    row.names = paste0("ENSG", seq_len(nrow(pbmc_small)))
  )

  vars <- Seurat::HVFInfo(object = pbmc_small, assay = "RNA", method = "vst")
  annotations <- pbmc_small@misc[["gene_annotations"]]
  vars$SYMBOL <- rownames(vars)
  vars$ENSEMBL <- annotations$input[match(rownames(vars), annotations$name)]
  pbmc_small@misc[["gene_dispersion"]] <- vars

  set.seed(0)
  pbmc_small[["percent.mt"]] <- rnorm(ncol(pbmc_small), 5, 1)
  pbmc_small[["doublet_scores"]] <- rnorm(ncol(pbmc_small), 0.5, 0.1)

  return(pbmc_small)
}

test_that("GetMitochondrialContent generates the expected return format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getMitochondrialContent(req, data)

  # number of values is number of cells
  expect_length(res, ncol(data))

  # result is a list of numeric values
  expect_type(res, "list")
  expect_type(unlist(res), "double")

  # result derived from percent.mt column
  expect_equal(unlist(res), unname(data$percent.mt[order(data$cells_id)]))
})

test_that("getDoubletScore generates the expected return format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getDoubletScore(req, data)

  # number of values is number of cells
  expect_length(res, ncol(data))

  # result is a list of numeric values
  expect_type(res, "list")
  expect_type(unlist(res), "double")

  # result derived from percent.mt column
  expect_equal(unlist(res), unname(data$doublet_scores[order(data$cells_id)]))
})

test_that("formatMetadataResult throws an error if metadata column is missing", {
  data <- mock_scdata()

  # these exist and are used
  expect_error(formatMetadataResult(data, column = "percent.mt"), NA)
  expect_error(formatMetadataResult(data, column = "doublet_scores"), NA)

  # but not this
  expect_error(formatMetadataResult(data, column = "blah"), "blah is not computed for this experiment.")
})

test_that("formatMetadataResult adds placeholders for filtered cells", {
  data <- mock_scdata()

  # remove 2 cells
  ncells.init <- ncol(data)
  exclude_idx <- c(5, 6)
  exclude_cell_ids <- data$cells_id[exclude_idx]
  data <- data[, -exclude_idx]
  expect_equal(ncol(data) + 2, ncells.init)

  res <- formatMetadataResult(data, column = "percent.mt")

  # result fills them in
  expect_length(res, ncells.init)

  # filled values are NULL
  for (cell_id in exclude_cell_ids){
    expect_null(res[[cell_id + 1]])
  }

  # rest are fine (not NULL)
  expect_type(unlist(res), "double")
  expect_length(unlist(res), ncol(data))
})


test_that("formatMetadataResult returns results in order of increasing cells_id", {
  data <- mock_scdata()

  # swap order of two cells
  ord <- new.ord <- seq_len(ncol(data))
  new.ord[c(3, 4)] <- ord[c(4, 3)]
  data@meta.data <- data@meta.data[new.ord, ]

  res <- formatMetadataResult(data, column = "percent.mt")

  # res not the same at percent.mt column (order is different)
  expect_false(identical(unlist(res), unname(data$percent.mt)))

  # if swap again will be the same
  data@meta.data <- data@meta.data[order(data$cells_id), ]
  expect_true(identical(unlist(res), unname(data$percent.mt)))
})


test_that("GetMitochondrialContent returns the same snapshot", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getMitochondrialContent(req, data)
  expect_snapshot(res)

})


test_that("getDoubletScore returns the same snapshot", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getDoubletScore(req, data)
  expect_snapshot(res)
})


test_that("getNGenes returns the same snapshot", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getNGenes(req, data)
  expect_snapshot(res)
})


test_that("getNUmis returns the same snapshot", {
  data <- mock_scdata()
  req <- mock_req()

  res <- getNUmis(req, data)
  expect_snapshot(res)
})
