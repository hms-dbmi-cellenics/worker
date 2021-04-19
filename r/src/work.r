library(RestRserve)
library(Matrix)
require(data.table)
library(RJSONIO)
library(Seurat)
library(sccore)

source("./differential_expression.r")
source("./embedding.r")
source("./get_metadata_information.r")
source("./expression.r")
source("./list_genes.r")
source("./cluster.r")

experiment_id <- Sys.getenv("EXPERIMENT_ID", unset = "e52b39624588791a7889e39c617f669e")

debug_step <- Sys.getenv("DEBUG_STEP", unset = "")
# over-ride manually to hot-reload
# debug_step <- "getClusters"

load_data <- function() {
    message(paste("Welcome to Biomage R worker, experiment id", experiment_id))

    loaded <- F
    data <- NULL

    while (!loaded) {
        data <- tryCatch(
            {
                print("Current working directory:")
                print(getwd())
                print("Experiment folder status:")
                print(list.files(paste("/data",experiment_id,sep = "/"),all.files=TRUE,full.names=TRUE))
                f <- readRDS(
                    paste(
                        "/data",experiment_id,"r.rds",
                        sep = "/"
                    )
                )
                loaded <- T
                length <- dim(f)
                message(
                    paste(
                        "Data successfully loaded, dimensions",
                        length[1], "x", length[2]
                    )
                )

                f
            },
            warning = function(w) {
                message(paste("file could not be loaded: ", w))
            },
            error = function(e) {
                message(paste("file could not be loaded: ", e))
            }
        )
        Sys.sleep(1)
    }

    return(data)
}

run_post <- function(req, post_fun) {
    handle_debug(req)
    post_fun(req)
}

handle_debug <- function(req) {
    task_name <- basename(req$path)
    is_debug <- debug_step == task_name
    
    if (is_debug) {
        fname <- paste0(task_name, '.RData')
        fpath_cont <- file.path('/debug', fname)
        fpath_host <- file.path('./data/debug', fname)
        
        message(sprintf('⚠️ DEBUG_STEP = %s. Saving `req` and `data` object.', task_name))
        # saving data for convenience
        # it (/data/e52../r.rds) is unchanged by worker
        save(data, req, file = fpath_cont)
        message(sprintf("⚠️ RUN load('%s') to restore environment.", fpath_host))
    }
}

create_app <- function(last_modified) {    
    last_modified_mw <- Middleware$new(
        process_request = function(request, response) {
            if (!file.info(path)$mtime == last_modified) {
                raise(
                    HTTPError$conflict(
                        body = toJSON(list(error = "The file is out of date and is currently being updated."))
                    )
                )
            }

            return(request)
        },
        id = "last_modified_mw"
    )

    app <- Application$new(
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
            result <- run_post(req, runDE)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getEmbedding",
        FUN = function(req, res) {
            result <- run_post(req, runEmbedding)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getDoubletScore",
        FUN = function(req, res) {
            result <- run_post(req, getDoubletScore)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getMitochondrialContent",
        FUN = function(req, res) {
            result <- run_post(req, getMitochondrialContent)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getExpression",
        FUN = function(req, res) {
            result <- run_post(req, runExpression)
            res$set_body(result)
    	}
    )
    app$add_post(
        path = "/v0/listGenes",
        FUN = function(req, res) {
            result <- run_post(req, getList)
            res$set_body(result)
    	}
    )
    app$add_post(
        path = "/v0/getClusters",
        FUN = function(req, res) {
            str(req$body)
            result <- run_post(req, getClusters)
            res$set_body(result)
    	}
    )

    return(app)
}

backend <- BackendRserve$new()

repeat {
    data <- load_data()
    path <- paste(
        "/data",experiment_id,"r.rds",
        sep = "/"
    )

    last_modified <- file.info(path)$mtime
    app <- create_app(last_modified)
    proc <- backend$start(app, http_port = 4000, background = TRUE)

    while(file.info(path)$mtime == last_modified) {
        Sys.sleep(10);
    }

    proc$kill()
}
