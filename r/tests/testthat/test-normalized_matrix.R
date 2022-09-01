mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  return(pbmc_small)
}

mock_req <- function(apply_filter = TRUE) {
  req <- list(
    body = list(
      applyFilter = apply_filter,
      filterBy = c(2:10, 30:36, 60:70)
        )
      )
}


test_that("GetNormalizedExpression generates the expected data frame format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- GetNormalizedExpression(req, data)
  expect_s3_class(res, "data.frame")
})

test_that("subsetting is applied and changes GetNormalizedExpression output", {
  data <- mock_scdata()
  req <- mock_req()

  res_filt <- GetNormalizedExpression(req, data)

  req <- mock_req(apply_filter = FALSE)
  res_unfilt <- GetNormalizedExpression(req, data)

  expect_false(identical(res_unfilt, res_filt))
})

test_that("GetNormalizedExpression correctly subsets the data", {
  data <- mock_scdata()
  req <- mock_req()

  subset_ids <- req$body$filterBy

  res <- GetNormalizedExpression(req, data)

  expect_false(ncol(data) == ncol(res))
  expect_equal(length(subset_ids), ncol(res))
})


test_that("GetNormalizedExpression doesn't subset the data when applyFilter is FALSE", {
  data <- mock_scdata()
  req <- mock_req(apply_filter = FALSE)

  res <- GetNormalizedExpression(req, data)

  expect_true(ncol(data) == ncol(res))
})
