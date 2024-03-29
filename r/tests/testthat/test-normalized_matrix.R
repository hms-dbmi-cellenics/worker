mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  return(pbmc_small)
}

mock_req <- function(apply_subset = TRUE) {
  subset_by <- c(2:10, 30:36, 60:70)

  req <- list(
    body = list(
      subsetBy = subset_by,
      applySubset = apply_subset
    )
  )
}


stub_fwrite <- function(matrix, fpath, quote = FALSE, row.names = TRUE) {
  # add '.' to path to point it to cwd instead of system root
  fpath <- paste0(".", fpath)
  if (!dir.exists(dirname(fpath))) {
    dir.create(dirname(fpath), recursive = TRUE)
  }
  data.table::fwrite(matrix, fpath, quote = FALSE, row.names = TRUE)
  return(fpath)
}



stubbed_GetNormalizedExpression <- function(req, data) {
  mockery::stub(
    GetNormalizedExpression,
    "data.table::fwrite",
    stub_fwrite
  )

  res <- GetNormalizedExpression(req, data)
  return(paste0(".", res))
}


test_that("GetNormalizedExpression saves the normalized matrix using the correct path", {
  data <- mock_scdata()
  req <- mock_req()

  res <- stubbed_GetNormalizedExpression(req, data)
  withr::defer(unlink(dirname(res), recursive = TRUE))

  expect_type(res, "character")
  expect_equal(gsub("^\\.", "", res), TMP_RESULTS_PATH_GZ)
})


test_that("GetNormalizedExpression correctly subsets the data", {
  data <- mock_scdata()
  req <- mock_req()

  original_ncol <- ncol(data)

  res <- stubbed_GetNormalizedExpression(req, data)
  withr::defer(unlink(dirname(res), recursive = TRUE))

  matrix <- data.table::fread(res, sep=",", quote = F, header = TRUE)

  expect_equal(ncol(matrix) - 2, length(req$body$subsetBy)) # -1 because of the rownames and index columns
  expect_false(ncol(data) == ncol(matrix))
})

test_that("subsetting is applied and changes GetNormalizedExpression output", {
  data <- mock_scdata()
  req <- mock_req()

  res_filt <- stubbed_GetNormalizedExpression(req, data)
  withr::defer(unlink(dirname(res_filt), recursive = TRUE))
  matrix_filt <- data.table::fread(res_filt, sep=",", quote = F, header = TRUE)

  req <- mock_req(apply_subset = FALSE)
  res_unfilt <- stubbed_GetNormalizedExpression(req, data)
  withr::defer(unlink(dirname(res_unfilt), recursive = TRUE))
  matrix_unfilt <- data.table::fread(res_unfilt, sep=",", quote = F, header = TRUE)

  expect_false(identical(matrix_unfilt, matrix_filt))
})


test_that("GetNormalizedExpression doesn't subset the data when applySubset is FALSE", {
  data <- mock_scdata()
  req <- mock_req(apply_subset = FALSE)

  expect_message(
    {
      res <- stubbed_GetNormalizedExpression(req, data)
      withr::defer(unlink(dirname(res), recursive = TRUE))
    },
    "No subsetting specified, sending the whole matrix"
  )

  matrix <-data.table::fread(res, sep=",", quote = F, header = TRUE)

  expect_true(ncol(data) == ncol(matrix) - 2)
})


test_that("GetNormalizedExpression fails if subsetBy is empty", {
  data <- mock_scdata()
  req <- mock_req()

  req$body$subsetBy <- c()

  expect_error(GetNormalizedExpression(req, data), "No cells match requested filters.")
})
