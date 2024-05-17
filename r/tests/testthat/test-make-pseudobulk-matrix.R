mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject", envir = environment())
  scdata <- pbmc_small

  scdata$samples <-
    rep(c("s1", "s2", "s3", "s4"), each = nrow(scdata@meta.data) / 4)

  # add annotations
  scdata@misc$gene_annotations <- data.frame(
    input = row.names(scdata),
    name = paste0("ENS", seq_len(nrow(scdata))),
    row.names = row.names(scdata)
  )

  # add custom sets
  sample <- scdata$samples
  scdata$custom <-
    dplyr::case_when(
      sample == "s1" ~ "base",
      sample == "s2" ~ "base",
      sample == "s3" ~ "background",
      sample == "s4" ~ "background"
    )

  # add cells not belonging to groups to compare
  scdata$custom[floor(ncol(scdata) / 10) * 1:10] <- NA

  return(scdata)
}


test_that("makePseudobulkMatrix returns object of correct type and dims", {
  scdata <- mock_scdata()
  res <- makePseudobulkMatrix(scdata)

  expected_rows <- nrow(scdata)
  expected_cols <-
    length(unique(scdata$samples[!is.na(scdata$custom)]))

  expect_equal(nrow(res), expected_rows)
  expect_equal(ncol(res), expected_cols)
  expect_s4_class(res, "Seurat")
})


test_that("makePseudobulkMatrix returns correct gene annotations", {
  scdata <- mock_scdata()
  res <- makePseudobulkMatrix(scdata)

  expected_gene_annot <- scdata@misc$gene_annotations

  expect_equal(res@misc$gene_annotations, expected_gene_annot)
})


test_that("makePseudobulkMatrix returns correct custom and samples slots", {
  scdata <- mock_scdata()
  res <- makePseudobulkMatrix(scdata)

  expected_metadata <- scdata@meta.data %>%
    dplyr::filter(!is.na(custom)) %>%
    dplyr::select(samples, custom) %>%
    unique() %>%
    `rownames<-`(.$samples)

  expect_equal(res@meta.data[, c("samples", "custom")], expected_metadata)
})


test_that("makePseudobulkMatrix correctly aggregates counts", {
  scdata <- mock_scdata()

  for (sample in unique(scdata$samples)) {
    sample_cells <- scdata@meta.data %>%
      dplyr::filter(!is.na(custom), samples == sample) %>%
      rownames()

    # aggregate manually,
    cell_counts <- scdata[["RNA"]]$counts[, sample_cells]
    expected_agg_counts <- Matrix::rowSums(cell_counts)

    res <- makePseudobulkMatrix(scdata)
    res_agg_counts <- res[["RNA"]]$counts[, sample]

    expect_equal(res_agg_counts, expected_agg_counts)
  }
})

test_that("makePseudobulkMatrix removes samples with fewer than 10 cells", {
  scdata <- mock_scdata()

  # doesn't filter samples with 10 or more cells
  res <- makePseudobulkMatrix(scdata)

  expected_cols <-
    length(unique(scdata$samples[!is.na(scdata$custom)]))
  expect_equal(ncol(res), expected_cols)


  # make s1 have 9 cells
  s1.cells <- colnames(scdata)[scdata$samples == 's1']
  filt.cells <- setdiff(s1.cells, head(s1.cells, 9))

  scdata_s1_9cells <- scdata[, !colnames(scdata) %in% filt.cells]
  expected_cols_s1_9cells <- sum(table(scdata_s1_9cells$samples) >= 10)
  expect_equal(expected_cols_s1_9cells, expected_cols-1)

  # s1 filtered
  res_s1_9cells <- makePseudobulkMatrix(scdata_s1_9cells)
  expect_equal(expected_cols_s1_9cells, ncol(res_s1_9cells))
})
