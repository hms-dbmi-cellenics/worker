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

create_app <- function(data) {
    app <- Application$new(content_type = "application/json")
    app$add_get(
        path = "/health",
        FUN = function(request, response) {
            response$set_body("up")
        }
    )
    app$add_post(
        path = "/v0/DifferentialExpression",
        FUN = function(req, res) {
            result <- runDE(req)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getEmbedding",
        FUN = function(req, res) {
            result <- runEmbedding(req)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getDoubletScore",
        FUN = function(req, res) {
            result <- getDoubletScore(req)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getMitochondrialContent",
        FUN = function(req, res) {
            result <- getMitochondrialContent(req)
            res$set_body(result)
        }
    )
    app$add_post(
        path = "/v0/getExpression",
        FUN = function(req, res) {
            result <- runExpression(req)
            res$set_body(result)
    	}
    )
    app$add_post(
        path = "/v0/listGenes",
        FUN = function(req, res) {
            result <- getList(req)
            res$set_body(result)
    	}
    )
    app$add_post(
        path = "/v0/getClusters",
        FUN = function(req, res) {
            result <- getClusters(req)
            res$set_body(result)
    	}
    )
    app$add_post(
        path = "/v0/reload",
        FUN = function(req, res) {
            system('/sbin/killall5')
            res$set_body('ok')
    	}
    )

    return(app)
}

repeat {
    data <- load_data()
    backend <- BackendRserve$new()
    app <- backend$start(create_app(data), http_port = 4000, background = TRUE)

    path <- paste(
        "/data",experiment_id,"r.rds",
        sep = "/"
    )
    
    last <- file.info(path)$mtime

    while(file.info(path)$mtime == last) {
        Sys.sleep(1);
    }

    app$kill()
}
