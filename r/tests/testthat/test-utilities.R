mock_scratchpad_cellset_object <- function(n) {
  # ensure cellIds is an int vector. same as when created by getExpressionCellSet
  list(
    key = "cell_set_key",
    name = "cell_set_name",
    rootNode = FALSE,
    color = "color",
    cellIds = as.integer(1:n)
  )
}

with_fake_http(test_that("sendCellSetToApi sends patch request", {
  scr_cellset_object <- mock_scratchpad_cellset_object(5)
  expect_PATCH(sendCellsetToApi(scr_cellset_object, "api_url", "experiment_id", "cell_set_key", "auth"))
}))



with_fake_http(test_that(
  "sendCellSetToApi sends the cellIds as an array, when there are >1 cells",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(10)
    req <- sendCellsetToApi(scr_cellset_object,
                            "api_url",
                            "experiment_id",
                            "cell_set_key",
                            "auth")

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    # extract the cellIds slot in the request
    cell_ids <-
      req_body[[1]]$`$match`$value$children[[1]]$`$insert`$value$cellIds

    expect_type(cell_ids, "list")

  }
))


with_fake_http(test_that(
  "sendCellSetToApi sends the cellIds as an array, when there is exactly 1 cell",
  {
    scr_cellset_object <- mock_scratchpad_cellset_object(1)
    req <- sendCellsetToApi(scr_cellset_object,
                            "api_url",
                            "experiment_id",
                            "cell_set_key",
                            "auth")

    # get req body, parsed as a json
    req_body <-
      httr::content(req, as = "parsed", type = "application/json")

    # extract the cellIds slot in the request
    cell_ids <-
      req_body[[1]]$`$match`$value$children[[1]]$`$insert`$value$cellIds

    expect_type(cell_ids, "list")

  }
))


test_that("complete_variable fills a vector with NAs where there are filtered cells", {

  mock_variable <- rnorm(50)
  mock_cell_ids <- seq.int(0,99, 2)

  res <- complete_variable(mock_variable, mock_cell_ids)

  expect_equal(length(res), max(mock_cell_ids) + 1)

})


test_that("complete_variable returns results ordered by increasing cell id", {
  set.seed(10)
  mock_variable <- rnorm(50)
  mock_cell_ids <- seq.int(0,99, 2)

  # shuffle cell ids
  mock_cell_ids <- sample(mock_cell_ids,size = length(mock_cell_ids), replace = FALSE)

  # sort the variable by cell id
  cell_id_sorted_variable <- mock_variable[order(mock_cell_ids)]

  res <- complete_variable(mock_variable, mock_cell_ids)

  # remove NAs from result vector
  expect_equal(res[!is.na(res)], cell_id_sorted_variable)

})
