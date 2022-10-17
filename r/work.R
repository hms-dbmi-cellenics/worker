library(Seurat)
library(dplyr)

for (f in list.files("R", ".R$", full.names = TRUE)) source(f)
load('R/sysdata.rda') # constants

load_data <- function(fpath) {
  loaded <- FALSE
  data <- NULL

  while (!loaded) {
    data <- tryCatch(
      {
        print("Current working directory:")
        print(getwd())
        print("Experiment folder status:")
        print(list.files(dirname(fpath), all.files = TRUE, full.names = TRUE))
        f <- readRDS(fpath)
        loaded <- TRUE
        length <- dim(f)
        message(
          "Data successfully loaded, dimensions",
          length[1], "x", length[2]
        )

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
    # over-ride manually to hot-reload
    # debug_step <- "getClusters"
    debug_step <- Sys.getenv("DEBUG_STEP", unset = "")

    handle_debug(req, debug_step)

    message(rep("✧",100))
    message("➥ Starting ",sub("run","",basename(req$path)))
    message("Input:")
    message(str(req$body))

    tryCatch({
            message("\nSeurat logs:")
            message("➡️ \n")
            tstart <- Sys.time()
            res <- post_fun(req, data)
            message("\n⬅️")

            message("\nResult length: ",length(res))
            message("\nResult head: ")
            message(str(head(res,10)))

            ttask <- format(Sys.time()-tstart, digits = 2)
            message("\n⏱️ Time to complete ", req$body$name, " for experiment ", experiment_id, ": ", ttask, '\n')
            message("✅ Finished ", req$body$name)
            message(rep("✧",100))

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

handle_debug <- function(req, debug_step) {
  task_name <- basename(req$path)
  is_debug <- debug_step == task_name | debug_step == "all"

  if (is_debug) {
    message(sprintf("⚠ DEBUG_STEP = %s. Saving `req` object.", task_name))
    req_fname <- sprintf("%s_%s_req.rds", experiment_id, task_name)
    saveRDS(req, file.path("/debug", req_fname))

    req_host <- file.path("./data/debug", req_fname)
    message(sprintf("⚠ RUN req  <- readRDS('%s') to restore 'req' object.", req_host))

    # copy data to /debug if doesn't exist
    data_fname <- sprintf("%s_data.rds", experiment_id)
    data_cont <- file.path("/debug", data_fname)

    if (!file.exists(data_cont)) {
      data_path <- file.path("/data", experiment_id, "r.rds")
      file.copy(data_path, data_cont)
    }

    data_host <- file.path("./data/debug", data_fname)
    message(sprintf("⚠ RUN data <- readRDS('%s') to restore 'data' object.", data_host))
  }
}

create_app <- function(last_modified, data, fpath) {
  last_modified_mw <- RestRserve::Middleware$new(
    process_request = function(request, response) {
      if (!file.info(fpath)$mtime == last_modified) {
        RestRserve::raise(
          RestRserve::HTTPError$conflict(
            body = RJSONIO::toJSON(list(error = "The file is out of date and is currently being updated."))
          )
        )
      }

      return(request)
    },
    id = "last_modified_mw"
  )

  app <- RestRserve::Application$new(
    content_type = "application/json"
  )

  app$append_middleware(last_modified_mw)

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
    path = "/v0/runExpression",
    FUN = function(req, res) {
      result <- run_post(req, runExpression, data)
      res$set_body(result)
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
    path = "/v0/runTrajectoryAnalysis",
    FUN = function(req, res) {
      result <- run_post(req, runTrajectoryAnalysis, data)
      res$set_body(result)
    }
  )
  app$add_post(
    path = "/v0/runGenerateTrajectoryGraph",
    FUN = function(req, res) {
      result <- run_post(req, runGenerateTrajectoryGraph, data)
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
fpath <- file.path("/data", experiment_id, "r.rds")

repeat {
  # need to load here as can change e.g. integration method
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
