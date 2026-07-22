# Tests for the CASSIA (Amazon Bedrock) annotation path. External calls -
# CASSIA's Python pipeline and the IAM bearer-token mint - are stubbed so the
# tests stay deterministic and offline.

# --- run_cassia -------------------------------------------------------------

# Minimal FindAllMarkers-style data frame with `n_clusters` clusters.
mock_markers <- function(n_clusters, genes_per = 3) {
  data.frame(
    cluster = rep(seq_len(n_clusters), each = genes_per),
    gene = paste0("GENE", seq_len(n_clusters * genes_per)),
    avg_log2FC = 1,
    stringsAsFactors = FALSE
  )
}

# Run run_cassia with the CASSIA pipeline + merge-back stubbed out, returning
# the arguments runCASSIA_pipeline was called with. The pipeline stub writes a
# fake results CSV into output_dir so run_cassia's downstream steps proceed.
capture_pipeline_args <- function(markers,
                                  provider = "https://bedrock/openai/v1",
                                  additional_info = NULL) {
  captured <- new.env()

  fake_pipeline <- function(...) {
    captured$args <- list(...)
    csv <- file.path(captured$args$output_dir, "results_FINAL_RESULTS.csv")
    utils::write.csv(data.frame(cluster = 1), csv, row.names = FALSE)
    invisible(NULL)
  }

  mockery::stub(run_cassia, "CASSIA::runCASSIA_pipeline", fake_pipeline)
  mockery::stub(
    run_cassia, "CASSIA::add_cassia_to_seurat",
    function(seurat_obj, ...) seurat_obj
  )

  run_cassia(
    data = "dummy_seurat",
    markers = markers,
    tissue = "Large Intestine",
    species = "Human",
    model = "qwen.qwen3-next-80b-a3b",
    provider = provider,
    additional_info = additional_info,
    cluster_col = "cassia_cluster"
  )

  captured$args
}

test_that("run_cassia runs one worker per cluster", {
  args <- capture_pipeline_args(mock_markers(4))
  expect_equal(args$max_workers, 4)
})

test_that("run_cassia caps max_workers at 50", {
  args <- capture_pipeline_args(mock_markers(60))
  expect_equal(args$max_workers, 50)
})

test_that("run_cassia forwards the provider, model and skips pre-flight validation", {
  args <- capture_pipeline_args(mock_markers(3), provider = "https://x/openai/v1")
  expect_equal(args$overall_provider, "https://x/openai/v1")
  expect_equal(args$annotation_model, "qwen.qwen3-next-80b-a3b")
  expect_equal(args$score_model, "qwen.qwen3-next-80b-a3b")
  expect_false(args$validate_api_keys_before_start)
})

test_that("run_cassia forwards non-empty additional_info", {
  args <- capture_pipeline_args(mock_markers(2), additional_info = "3 tumor, 2 normal")
  expect_equal(args$additional_info, "3 tumor, 2 normal")
})

test_that("run_cassia treats blank additional_info as no context", {
  args <- capture_pipeline_args(mock_markers(2), additional_info = "   ")
  expect_null(args$additional_info)
})

# --- set_cassia_bedrock_auth ------------------------------------------------

# Fake reticulate::import returning the token generator and an `os` whose
# environ$update records what was written. `recorder` captures both.
fake_reticulate_import <- function(recorder) {
  function(module) {
    if (module == "aws_bedrock_token_generator") {
      list(provide_token = function(region) {
        recorder$region <- region
        "bedrock-api-key-FAKE"
      })
    } else if (module == "os") {
      list(environ = list(update = function(x) recorder$env <- x))
    } else {
      stop("unexpected module: ", module)
    }
  }
}

test_that("set_cassia_bedrock_auth mints a token into CUSTOMIZED_API_KEY and returns the endpoint", {
  withr::local_envvar(AWS_REGION = "", AWS_DEFAULT_REGION = "")
  rec <- new.env()
  mockery::stub(
    set_cassia_bedrock_auth, "reticulate::import", fake_reticulate_import(rec)
  )

  url <- set_cassia_bedrock_auth()

  expect_equal(rec$region, "us-east-1")
  expect_equal(rec$env$CUSTOMIZED_API_KEY, "bedrock-api-key-FAKE")
  expect_equal(url, "https://bedrock-runtime.us-east-1.amazonaws.com/openai/v1")
})

test_that("set_cassia_bedrock_auth honours AWS_REGION", {
  withr::local_envvar(AWS_REGION = "eu-west-1")
  rec <- new.env()
  mockery::stub(
    set_cassia_bedrock_auth, "reticulate::import", fake_reticulate_import(rec)
  )

  url <- set_cassia_bedrock_auth()

  expect_equal(rec$region, "eu-west-1")
  expect_match(url, "bedrock-runtime\\.eu-west-1\\.amazonaws\\.com")
})

test_that("set_cassia_bedrock_auth raises a CASSIA error when the token can't be minted", {
  failing_import <- function(module) {
    if (module == "aws_bedrock_token_generator") stop("no credentials")
    list(environ = list(update = function(x) NULL))
  }
  mockery::stub(set_cassia_bedrock_auth, "reticulate::import", failing_import)

  expect_error(set_cassia_bedrock_auth(), "Bedrock")
})
