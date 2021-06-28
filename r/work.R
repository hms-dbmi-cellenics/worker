library(Seurat)
library(dplyr)

for (f in list.files('R', '.R$', full.names = TRUE)) source(f)

load_data <- function(fpath) {

    loaded <- FALSE
    data <- NULL

    while (!loaded) {
        data <- tryCatch(
            {
                print("Current working directory:")
                print(getwd())
                print("Experiment folder status:")
                print(list.files(dirname(fpath), all.files=TRUE, full.names=TRUE))
                f <- readRDS(fpath)
                loaded <- TRUE
                length <- dim(f)
                message("Data successfully loaded, dimensions",
                        length[1], "x", length[2])

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
    post_fun(req, data)
}

handle_debug <- function(req, debug_step) {
    task_name <- basename(req$path)
    is_debug <- debug_step == task_name | debug_step == 'all'

    if (is_debug) {
        message(sprintf('⚠ DEBUG_STEP = %s. Saving `req` object.', task_name))
        req_fname <- sprintf('%s_%s_req.rds', experiment_id, task_name)
        saveRDS(req, file.path('/debug', req_fname))

        req_host  <- file.path('./data/debug', req_fname)
        message(sprintf("⚠ RUN req  <- readRDS('%s') to restore 'req' object.",  req_host))

        # copy data to /debug if doesn't exist
        data_fname <- sprintf('%s_data.rds', experiment_id)
        data_cont <- file.path('/debug', data_fname)

        if (!file.exists(data_cont)) {
            data_path <- file.path('/data', experiment_id, 'r.rds')
            file.copy(data_path, data_cont)
        }

        data_host <- file.path('./data/debug', data_fname)
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
        path = "/v0/getExpression",
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
            str(req$body)
            result <- run_post(req, getClusters, data)
            res$set_body(result)
        }
    )

    return(app)
}

experiment_id <- Sys.getenv("EXPERIMENT_ID", unset = "e52b39624588791a7889e39c617f669e")
message(paste("Welcome to Biomage R worker, experiment id", experiment_id))

backend <- RestRserve::BackendRserve$new()
fpath <- file.path("/data", experiment_id, "r.rds")

repeat {
    # need to load here as can change e.g. integration method
    data <- load_data(fpath)
    last_modified <- file.info(fpath)$mtime
    app <- create_app(last_modified, data, fpath)
    proc <- backend$start(app, http_port = 4000, background = TRUE)

    while(file.info(fpath)$mtime == last_modified) {
        Sys.sleep(10);
    }

    proc$kill()
}

