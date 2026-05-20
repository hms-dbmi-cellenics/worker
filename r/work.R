
suppressPackageStartupMessages( {
  library(Seurat)
  library(dplyr)
})

# v5 is default but making explicit
options(Seurat.object.assay.version = "v5")

for (f in list.files("R", ".R$", full.names = TRUE)) source(f)
load("R/sysdata.rda") # constants

read_flex <- function(fpath) {
  data_dir <- dirname(fpath)

  if (tools::file_ext(fpath) == "qs") {
    data <- qs::qread(fpath)
  } else {
    data <- readRDS(fpath)
  }

  is_bpcells <- is(data[["RNA"]]$counts, "IterableMatrix")
  if (is_bpcells) {
    data <- load_bpcells(data, data_dir)
  }

  return(data)
}

load_data <- function(fpath) {
  loaded <- FALSE
  data <- NULL

  while (!loaded) {
    data <- tryCatch(
      {
        message(
          "\nCurrent working directory: ", getwd(),
          "\nExperiment folder status:\n",
          paste(
            list.files(dirname(fpath), all.files = TRUE, full.names = TRUE),
            collapse = "\n"
          )
        )

        f <- read_flex(fpath)
        loaded <- TRUE
        length <- dim(f)

        message(
          "\nData successfully loaded, dimensions: ",
          length[1], " x ", length[2]
        )

        message("\nSession info:")
        print(sessionInfo(),  locale = FALSE)
        cat("\n")

        return(f)
      },
      warning = function(w) {
        message("file could not be loaded: ", w)
      },
      error = function(e) {
        message("file could not be loaded: ", e)
      }
    )
    Sys.sleep(1)
  }

  return(data)
}

run_post <- function(req, post_fun, data) {

  message(rep("✧",100))
  message("➥ Starting ",sub("run","",basename(req$path)))
  message("Input:")
  str(req$body, list.len = 6)
  cat("\n")

  tryCatch({
    tstart <- Sys.time()
    res <- post_fun(req, data)

    message("\nResult length: ",length(res))
    message("Result head: ")
    str(head(res, 10))

    ttask <- format(Sys.time() - tstart, digits = 2)
    message(
      "\n⏱️ Time to complete ", req$body$name,
      " for experiment ", experiment_id, ": ", ttask, "\n"
    )
    message("✅ Finished ", req$body$name)
    message(rep("✧", 100))

    return(
      formatResponse(
        res,
        NULL
      )
    )
  },
  error = function(e) {
    message("🚩 --------- 🚩")
    message("Error at worker task: ", e$message)

    return(
      formatResponse(
        NULL,
        extractErrorList(e$message)
      )
    )
  }
  )
}

create_app <- function(last_modified, data, fpath) {
  last_modified_mw <- RestRserve::Middleware$new(
    process_request = function(request, response) {
      if (file.info(fpath)$mtime != last_modified) {
        RestRserve::raise(
          RestRserve::HTTPError$conflict(
            body = RJSONIO::toJSON(
              list(
                error = paste0(
                  "The file is out of date and is currently ",
                  "being updated."
                )
              )
            )
          )
        )
      }

      return(request)
    },
    id = "last_modified_mw"
  )

  encode_decode_middleware <- RestRserve::EncodeDecodeMiddleware$new()

  # the json encoder by default is not precise enough
  # so we set a custom one without precision limit (digits=NA)
  # yyjsonr is 2-10x faster than jsonlite
  encode_decode_middleware$ContentHandlers$set_encode(
    "application/json",
    function(x, unbox = TRUE)  {
      yyjsonr::write_json_str(
        x,
        opts = list(
          dataframe = "columns",
          digits_signif = 4,
          auto_unbox = unbox
        )
      )
    }
  )

  app <- RestRserve::Application$new(
    content_type = "application/json",
    middleware = list(encode_decode_middleware, last_modified_mw)
  )

  app$add_get(
    path = "/health",
    FUN = function(request, response) {
      response$set_body("up")
    }
  )
  app$add_post(
    path = "/v0/DifferentialExpression",
    FUN = function(req, res) {
      result <- run_post(req, runDE, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getEmbedding",
    FUN = function(req, res) {
      result <- run_post(req, runEmbedding, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getDoubletScore",
    FUN = function(req, res) {
      result <- run_post(req, getDoubletScore, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getMitochondrialContent",
    FUN = function(req, res) {
      result <- run_post(req, getMitochondrialContent, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getNGenes",
    FUN = function(req, res) {
      result <- run_post(req, getNGenes, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getNUmis",
    FUN = function(req, res) {
      result <- run_post(req, getNUmis, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runExpression",
    FUN = function(req, res) {
      result <- run_post(req, runExpression, data)
      res$set_body(result$data)
    }
  )
  app$add_post(
    path = "/v0/listGenes",
    FUN = function(req, res) {
      result <- run_post(req, getList, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getClusters",
    FUN = function(req, res) {
      result <- run_post(req, runClusters, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runMarkerHeatmap",
    FUN = function(req, res) {
      result <- run_post(req, runMarkerHeatmap, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getBackgroundExpressedGenes",
    FUN = function(req, res) {
      result <- run_post(req, getBackgroundExpressedGenes, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/getExpressionCellSet",
    FUN = function(req, res) {
      result <- run_post(req, getExpressionCellSet, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runTrajectoryAnalysisPseudoTimeTask",
    FUN = function(req, res) {
      result <- run_post(req, runTrajectoryAnalysisPseudoTimeTask, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runTrajectoryAnalysisStartingNodesTask",
    FUN = function(req, res) {
      result <- run_post(req, runTrajectoryAnalysisStartingNodesTask, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runDotPlot",
    FUN = function(req, res) {
      result <- run_post(req, runDotPlot, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/GetNormalizedExpression",
    FUN = function(req, res) {
      result <- run_post(req, GetNormalizedExpression, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/ScTypeAnnotate",
    FUN = function(req, res) {
      result <- run_post(req, ScTypeAnnotate, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/DownloadAnnotSeuratObject",
    FUN = function(req, res) {
      result <- run_post(req, DownloadAnnotSeuratObject, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/CellCycleScoring",
    FUN = function(req, res) {
      result <- run_post(req, cellCycleScoring, data)
      res$set_body(result)
    }
  )
  return(app)
}

repeat {
  label_path <- "/etc/podinfo/labels"
  experiment_id <- NA

  if (file.exists(label_path)) {
    labels <- read.csv(label_path, sep = "=", row.names = 1, header = FALSE)
    experiment_id <- labels["experimentId", ]
  }

  if (is.na(experiment_id)) {
    experiment_id <- Sys.getenv("EXPERIMENT_ID", unset = NA)
  }

  if (is.na(experiment_id)) {
    message("No experiment ID label set yet, waiting...")
    Sys.sleep(5)
  } else {
    message(paste("Welcome to Cellenics R worker, experiment id", experiment_id))
    break
  }
}

backend <- RestRserve::BackendRserve$new()

get_fpath <- function(experiment_id) {
  data_dir <- file.path("/data", experiment_id)

  fpath <- NULL

  while (is.null(fpath)) {
    fnames <- list.files(data_dir)

    if ("r.qs" %in% fnames) {
      fpath <- file.path(data_dir, "r.qs")

    } else if ("r.rds" %in% fnames) {
      # fallback to rds for older experiments
      fpath <- file.path(data_dir, "r.rds")
    }
    Sys.sleep(1)
  }

  return(fpath)
}

repeat {
  # need to load here as can change e.g. integration method
  cleanupMarkersCache()

  fpath <- get_fpath(experiment_id)

  data <- load_data(fpath)
  last_modified <- file.info(fpath)$mtime
  app <- create_app(last_modified, data, fpath)
  proc <- backend$start(app, http_port = 4000, background = TRUE)

  while (file.info(fpath)$mtime == last_modified) {
    Sys.sleep(10)
  }
  message("Detected a change in the rds object, reloading...")
  proc$kill()
}
