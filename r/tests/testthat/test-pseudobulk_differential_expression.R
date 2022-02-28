mock_pbulk <- function() {

  # mock pseudobulk SeuratObject with annotation
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )


  pbmc_raw <- Matrix::Matrix(as.matrix(pbmc_raw), sparse = TRUE)
  row.names(pbmc_raw) <- paste0("ENSG", seq_len(nrow(pbmc_raw)))

  groups <- rep(letters[1:4], each = 20)
  pbulk <- Matrix.utils::aggregate.Matrix(Matrix::t(pbmc_raw), groupings = groups, fun = "sum")
  pbulk <- Matrix::t(pbulk)

  return(pbulk)
}

# separated so that can add fake up/down-regulated genes
pbulk_to_seurat <- function(pbulk) {
  pbulk <- SeuratObject::CreateSeuratObject(counts = pbulk)
  pbulk@misc$gene_annotations <- data.frame(
    input = row.names(pbulk),
    name = paste0("SYMBOL", seq_len(nrow(pbulk))),
    row.names = row.names(pbulk)
  )

  # add comparison group
  pbulk$custom <- c("base", "base", "background", "background")


  return(pbulk)
}

test_that("runPseudobulkDE runs", {
  pbulk <- mock_pbulk()
  pbulk <- pbulk_to_seurat(pbulk)

  res <- runPseudobulkDE(pbulk)
  expect_equal(class(res), "data.frame")
})


test_that("the result of runPseudobulkDE detects up-regulated and down-regulated genes", {
  pbulk <- mock_pbulk()

  # add up-regulated and down-regulated gene
  fake <- matrix(c(
    564, 562, 24, 10,
    2, 5, 121, 130
  ),
  byrow = TRUE, nrow = 2, dimnames = list(c("UP", "DOWN"))
  )

  pbulk <- rbind(fake, pbulk)
  pbulk <- pbulk_to_seurat(pbulk)
  res <- runPseudobulkDE(pbulk)

  # correct direction
  expect_true(res["UP", "logFC"] > 0)
  expect_true(res["DOWN", "logFC"] < 0)

  # highest significance
  expect_setequal(head(row.names(res), 2), c("UP", "DOWN"))
})



test_that("the result of runPseudobulkDE includes all genes, including filtered ones", {
  pbulk <- mock_pbulk()
  pbulk <- pbulk_to_seurat(pbulk)

  res <- runPseudobulkDE(pbulk)

  # all ensemble ids result
  expect_setequal(row.names(res), row.names(pbulk))

  # some values in the table are NA (filtered genes)
  expect_true(anyNA(res))
})

test_that("runPseudobulkDE works for 1 sample vs 1 sample comparisons", {
  pbulk <- mock_pbulk()
  pbulk <- pbulk_to_seurat(pbulk)

  # make 1 vs 1 sample
  pbulk <- pbulk[, c(1, 3)]
  expect_setequal(pbulk$custom, c("base", "background"))

  res <- runPseudobulkDE(pbulk)

  # all ensemble ids in result
  expect_setequal(row.names(res), pbulk@misc$gene_annotations$input)

  # all expected columns in result
  expect_setequal(colnames(res), c("logFC", "AveExpr"))
})

test_that("runPseudobulkDE sorts by absolute logFC for 1 sample vs 1 sample comparisons", {
  pbulk <- mock_pbulk()

  # add up-regulated and down-regulated gene
  fake <- matrix(c(
    564, 562, 1, 1,
    2, 2, 1210, 1300
  ),
  byrow = TRUE, nrow = 2, dimnames = list(c("UP", "DOWN"))
  )

  pbulk <- rbind(fake, pbulk)
  pbulk <- pbulk_to_seurat(pbulk)

  # make 1 vs 1 sample
  pbulk <- pbulk[, c(1, 3)]
  expect_setequal(pbulk$custom, c("base", "background"))

  res <- runPseudobulkDE(pbulk)

  # correct direction
  expect_true(res["UP", "logFC"] > 0)
  expect_true(res["DOWN", "logFC"] < 0)

  # sorted by logFC
  expect_equal(res$name, res$name[order(abs(res$logFC), decreasing = TRUE)])
})
