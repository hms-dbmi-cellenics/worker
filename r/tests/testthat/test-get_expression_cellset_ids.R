test_that("getExpressionCellSetIDs with a single expression filter returns correct cell ids", {
  data <- mock_scdata()
  req <- list(list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5))
  res <- getExpressionCellSetIDs(req, data)

  annot <- data@misc$gene_annotations
  enid <- annot$input[annot$name == req[[1]]$geneName]
  expected <- data$cells_id[data[["RNA"]]$data[enid, ] > 0.5]
  expect_equal(res$keep_ids, expected)
})

test_that("getExpressionCellSetIDs runs with bpcells", {
  data <- mock_scdata(use_bpcells = TRUE)
  req <- list(list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5))
  res <- getExpressionCellSetIDs(req, data)

  annot <- data@misc$gene_annotations
  enid <- annot$input[annot$name == req[[1]]$geneName]
  data_mat <- as(data[["RNA"]]$data, "dgCMatrix")
  expected <- data$cells_id[data_mat[enid, ] > 0.5]
  expect_equal(res$keep_ids, expected)
})

test_that("getExpressionCellSetIDs with two expression filters returns correct cell ids", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0.5)
  )
  res <- getExpressionCellSetIDs(req, data)

  # gives same results as manual
  annot <- data@misc$gene_annotations
  enid1 <- annot$input[annot$name == req[[1]]$geneName]
  enid2 <- annot$input[annot$name == req[[2]]$geneName]
  expected <- data$cells_id[data[["RNA"]]$data[enid1, ] > 0.5 & data[["RNA"]]$data[enid2, ] > 0.5]

  expect_equal(res$keep_ids, expected)

  # order invariant
  req2 <- list(
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5)
  )

  res2 <- getExpressionCellSetIDs(req2, data)
  expect_equal(res$keep_ids, res2$keep_ids)
})


test_that("getExpressionCellSetIDs returns no cell ids if filters remove everything", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "MS4A1", comparisonType = "lessThan", thresholdValue = 0.5)
  )
  res <- getExpressionCellSetIDs(req, data)
  expect_length(res$keep_ids, 0)
})


test_that("getExpressionCellSetIDs returns all cell ids if filters include everything", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = -1)
  )
  res <- getExpressionCellSetIDs(req, data)
  expect_equal(res$keep_ids, data$cells_id)
})



test_that("getExpressionCellSetIDs fails if requested geneNames are not present in SeuratObject", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "BLAH", comparisonType = "greaterThan", thresholdValue = 0)
  )

  expect_error(getExpressionCellSetIDs(req, data), "gene name\\(s\\) that are not present[.]")
})

test_that("getExpressionCellSet fails if filters remove all cells", {
  data <- mock_scdata()
  genes_config <- list(
    list(geneName = "MS4A1", comparisonType = "greaterThan", thresholdValue = 0.5),
    list(geneName = "MS4A1", comparisonType = "lessThan", thresholdValue = 0.5)
  )

  req <- list(body = list(genesConfig = genes_config))
  expect_error(getExpressionCellSet(req, data), "No cells match requested filters.")
})


test_that("CellSet naming is correct", {
  data <- mock_scdata()
  req <- list(
    list(geneName = "CD79B", comparisonType = "greaterThan", thresholdValue = 0),
    list(geneName = "MS4A1", comparisonType = "lessThan", thresholdValue = 0.5)
  )

  res <- getExpressionCellSetIDs(req, data)

  expect_equal(res$cell_set_name, "CD79B>0, MS4A1<0.5")
})

test_that("We attempt to patch the API.", {
  data <- mock_scdata()
  req <- list(
    body = list(
      genesConfig = list(
        list(
          geneName = "MS4A1",
          comparisonType = "greaterThan",
          thresholdValue = 0.5
        )
      ),
      config = list(
        apiUrl = "http://host.docker.internal:3000",
        experimentId = "12345",
        authJwt = "1234"
      )
    )
  )

  expect_error(
    getExpressionCellSet(req, data),
    "host.docker.internal"
  )
})


with_fake_http(test_that("getExpressionCellSet sends patch request", {
  data <- mock_scdata()
  req <- list(
    body = list(
      genesConfig = list(
        list(
          geneName = "MS4A1",
          comparisonType = "greaterThan",
          thresholdValue = 0.5
        )
      ),
      config = list(
        apiUrl = "http://host.docker.internal:3000",
        experimentId = "12345",
        authJwt = "1234"
      )
    )
  )
  expect_PATCH(getExpressionCellSet(req, data))
}))
