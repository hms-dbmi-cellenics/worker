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

load_data <- function() {
    experiment_id <- Sys.getenv("EXPERIMENT_ID", unset = "5928a56c7cbff9de78974ab50765ed20")
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

    return(app)
}

data <- load_data()
backend <- BackendRserve$new()

backend$start(create_app(data), http_port = 4000)
