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

string_to_df <- function(text_df) {
  df <- read.table(text = text_df, sep =",", header=TRUE)

  # This part to rollback: rownames_to_column
  rownames(df) <- df[,1]
  df[,1] <- NULL

  return(df)
}

# test_that("GetNormalizedExpression generates the expected string format", {
#   data <- mock_scdata()
#   req <- mock_req()
#
#   res <- GetNormalizedExpression(req, data)
#   expect_type(res, "character")
# })

test_that("subsetting is applied and changes GetNormalizedExpression output", {
  data <- mock_scdata()
  req <- mock_req()

  res_filt <- GetNormalizedExpression(req, data)

  req <- mock_req(apply_subset = FALSE)
  res_unfilt <- GetNormalizedExpression(req, data)

  expect_false(identical(res_unfilt, res_filt))
})

# test_that("GetNormalizedExpression correctly subsets the data", {
#   data <- mock_scdata()
#   req <- mock_req()
#
#   subset_ids <- req$body$subsetBy
#
#   res <- GetNormalizedExpression(req, data)
#
#   df <- string_to_df(res)
#
#   expect_false(ncol(data) == ncol(df))
#   expect_equal(ncol(df), length(subset_ids))
# })


# test_that("GetNormalizedExpression doesn't subset the data when applySubset is FALSE", {
#   data <- mock_scdata()
#   req <- mock_req(apply_subset = FALSE)
#
#   res <- GetNormalizedExpression(req, data)
#
#   df <- string_to_df(res)
#
#   expect_true(ncol(data) == ncol(df))
# })


test_that("GetNormalizedExpression fails if subsetBy is empty", {
  data <- mock_scdata()
  req <- mock_req()

  req$body$subsetBy <- c()

  expect_error(GetNormalizedExpression(req, data), "No cells match requested filters.")
})

