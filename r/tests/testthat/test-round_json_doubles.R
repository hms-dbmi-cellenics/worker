test_that("round_json_doubles keeps 6 significant digits of non-whole doubles", {
  expect_equal(round_json_doubles(c(-3.456789123, 12345.6789, 1.23456789e-300)),
               c(-3.45679, 12345.7, 1.23457e-300))
})

test_that("round_json_doubles leaves whole-number doubles (e.g. cell ids) exact", {
  ids <- c(0, 123456, 788123, 1234567, 98765432)
  expect_identical(round_json_doubles(ids), ids)
  expect_identical(round_json_doubles(c(1.5, 1234567)), c(1.5, 1234567))
})

test_that("round_json_doubles leaves non-finite values, integers and other types alone", {
  expect_identical(round_json_doubles(c(NA, NaN, Inf, -Inf, 0.123456789)),
                   c(NA, NaN, Inf, -Inf, 0.123457))
  expect_identical(round_json_doubles(c(1L, 123456789L)), c(1L, 123456789L))
  expect_identical(round_json_doubles("0.123456789"), "0.123456789")
  expect_identical(round_json_doubles(TRUE), TRUE)
  expect_identical(round_json_doubles(NULL), NULL)
  date <- as.Date("2026-10-17")
  expect_identical(round_json_doubles(date), date)
})

test_that("round_json_doubles recurses into lists, data.frames and matrices keeping structure", {
  result <- list(
    cellIds = c(10, 20),
    embedding = list(NULL, c(1.23456789, 2.5), NULL),
    stats = data.frame(gene = c("A", "B"), logFC = c(0.123456789, 2), count = c(3L, 4L)),
    mat = matrix(c(1.23456789, 2, 3.33333333, 4), nrow = 2),
    nested = list(list(p = 0.000123456789))
  )

  out <- round_json_doubles(result)

  expect_identical(out$cellIds, c(10, 20))
  expect_null(out$embedding[[1]])
  expect_equal(out$embedding[[2]], c(1.23457, 2.5))
  expect_length(out$embedding, 3)
  expect_s3_class(out$stats, "data.frame")
  expect_equal(out$stats$logFC, c(0.123457, 2))
  expect_identical(out$stats$count, c(3L, 4L))
  expect_identical(out$stats$gene, c("A", "B"))
  expect_equal(dim(out$mat), c(2L, 2L))
  expect_equal(out$mat[, 2], c(3.33333, 4))
  expect_equal(out$nested[[1]]$p, 0.000123457)
})

test_that("round_json_doubles output encodes to shorter JSON with yyjsonr", {
  x <- list(c(-3.456789123456, 12.3456789123), NULL, c(1234, 5678))
  json <- yyjsonr::write_json_str(round_json_doubles(x), opts = list(auto_unbox = TRUE))
  expect_equal(json, "[[-3.45679,12.3457],[],[1234.0,5678.0]]")
})
