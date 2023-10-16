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

mockery::stub(
  GetNormalizedExpression,
  "vroom::vroom_write",
  INTERNAL_RESULTS_PATH
)

stubbed_GetNormalizedExpression <- function(req, data) {
  GetNormalizedExpression(req, data)
}

mock_GetNormalizedExpression <- function(req, data) {
  subset_ids <- req$body$subsetBy
  apply_subset <- req$body$applySubset

  if (apply_subset && length(subset_ids) == 0) {
    stop(
      generateErrorMessage(
        error_codes$EMPTY_CELL_SET,
        "No cells match requested filters."
      )
    )
  }

  message("Extracting normalized expression matrix")

  if (apply_subset) {
    message("Number of cells before subsetting: ", ncol(data))
    data <- subsetIds(data, cells_id = subset_ids)
  } else {
    message("No subsetting specified, sending the whole matrix")
  }

  matrix <- as.data.frame(Seurat::GetAssayData(data, slot = "data", assay = "RNA"))

  message("Number of cells in matrix to return: ", ncol(matrix))

  matrix <- tibble::rownames_to_column(matrix, var = " ")

  temp_path <- tempfile(fileext = ".csv")
  vroom::vroom_write(matrix, temp_path, delim = ",", quote = "none")
  return(temp_path)
}



test_that("GetNormalizedExpression saves the normalized matrix using the correct path", {
  data <- mock_scdata()
  req <- mock_req()

  res <- stubbed_GetNormalizedExpression(req, data)

  expect_type(res, "character")
  expect_equal(res, INTERNAL_RESULTS_PATH)
})


test_that("GetNormalizedExpression correctly subsets the data", {
  data <- mock_scdata()
  req <- mock_req()

  original_ncol <- ncol(data)

  res <- mock_GetNormalizedExpression(req, data)

  matrix <- vroom::vroom(res, delim = ",", col_names = TRUE)

  expect_equal(ncol(matrix) - 1, length(req$body$subsetBy)) # -1 because of the rownames column
  expect_false(ncol(data) == ncol(matrix))

  unlink(res)
})

test_that("subsetting is applied and changes GetNormalizedExpression output", {
  data <- mock_scdata()
  req <- mock_req()

  res_filt <- mock_GetNormalizedExpression(req, data)
  matrix_filt <- vroom::vroom(res_filt, delim = ",", col_names = TRUE)

  req <- mock_req(apply_subset = FALSE)
  res_unfilt <- mock_GetNormalizedExpression(req, data)
  matrix_unfilt <- vroom::vroom(res_unfilt, delim = ",", col_names = TRUE)

  expect_false(identical(matrix_unfilt, matrix_filt))

  unlink(res_filt)
  unlink(res_unfilt)
})


test_that("GetNormalizedExpression doesn't subset the data when applySubset is FALSE", {
  data <- mock_scdata()
  req <- mock_req(apply_subset = FALSE)

  expect_message({
    res <- mock_GetNormalizedExpression(req, data)
  }, "No subsetting specified, sending the whole matrix")

  matrix <- vroom::vroom(res, delim = ",", col_names = TRUE)

  expect_true(ncol(data) == ncol(matrix) - 1)

  unlink(res)
})


test_that("GetNormalizedExpression fails if subsetBy is empty", {
  data <- mock_scdata()
  req <- mock_req()

  req$body$subsetBy <- c()

  expect_error(GetNormalizedExpression(req, data), "No cells match requested filters.")
})

