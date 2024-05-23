mock_req <- function() {
  req <- list(
    body = list(
      useMarkerGenes = TRUE,
      customGenesList = c("PF4", "CST3"),
      numberOfMarkers = 5,
      applyFilter = FALSE,
      filterBy = list(
        cellIds = c(0:25, 51:75)
      ),
      groupBy = list(
        children = list(
          list(
            name = "hello",
            cellIds = 0:50
          ),
          list(
            name = "hello2",
            cellIds = 51:100
          )
        )
      )
    )
  )
}

mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = row.names(pbmc_small),
    name = row.names(pbmc_small),
    row.names = row.names(pbmc_small)
  )
  return(pbmc_small)
}

test_that("dotplot generates the expected list format", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDotPlot(req, data)
  expect_named(res, c("cellSetsIdx", "geneNameIdx", "avgExpression", "cellsPercentage", "cellSetsNames", "geneNames"))
  expect_type(res$cellSetsNames, "character")
  expect_type(res$geneNames, "character")
  expect_type(res$cellSetsIdx, "double")
  expect_type(res$geneNameIdx, "double")
  expect_type(res$avgExpression, "double")
  expect_type(res$cellsPercentage, "double")
})

test_that("useMarkerGenes works", {
  data <- mock_scdata()
  req <- mock_req()

  res <- runDotPlot(req, data)
  expect_snapshot(res)
})

test_that("customGenesList is used if useMarkerGenes is FALSE", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$useMarkerGenes <- FALSE

  res <- runDotPlot(req, data)
  genes_used <- res$geneNames

  expect_true(all(genes_used %in% req$body$customGenesList))
})

test_that("subsetting is applied and changes dotpot output", {
  data <- mock_scdata()
  req <- mock_req()

  res_unfilt <- runDotPlot(req, data)

  req$body$applyFilter <- TRUE
  res_filt <- runDotPlot(req, data)

  expect_false(identical(res_unfilt, res_filt))
})

test_that("Dotplot returns the correct values", {
  data <- mock_scdata()
  req <- mock_req()
  req$body$useMarkerGenes <- FALSE

  dotPlot_res <- runDotPlot(req, data)
  correct_res <- data.frame()

  for (group in req$body$groupBy$children) {
    group_counts <- data@assays$RNA$counts[req$body$customGenesList, intersect(group$cellIds, data$cells_id)]
    for (gene in rownames(group_counts)) {
      correct_res[group$name, gene] <- mean(expm1(group_counts[gene, ]))
    }
  }

  scaled_res <- scale(correct_res)
  scaled_res <- Seurat::MinMax(scaled_res, min = -2.5, max = 2.5)
  scaled_res <- unlist(as.list(x = t(x = scaled_res)))
  expect_equal(dotPlot_res$avgExpression, scaled_res)
})
