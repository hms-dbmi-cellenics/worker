library(RestRserve)
library(pagoda2)
library(conos)
library(Matrix)
require(data.table)
library(DESeq2)
library(RJSONIO)

experiment_id <- Sys.getenv("EXPERIMENT_ID")
message(paste("Welcome to Biomage R worker, experiment id", experiment_id))

data <- readRDS(paste("/data", experiment_id, "r.rds", sep = "/"))
length <- dim(data$counts)

message(
    paste("Data successfully loaded, dimensions", length[1], "x", length[2])
)

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

        # set up a factor with the appropriate cells
        factor <- rep(
            c("base", "background"),
            c(
                length(req$body$baseCells),
                length(req$body$backgroundCells)
            )
        )

        names(factor) <- c(
            req$body$baseCells, req$body$backgroundCells
        )

        data$clusters$custom$de <- factor(factor)

        # # compute differential expression
        data$getDifferentialGenes(
            type = "custom", verbose = T, clusterType = "de"
        )

        # # get the necessary results
        result <- data$diffgenes$custom$de$base

        # replace name with Gene names
        result$Gene <- data$misc$gene_annotations[
            match(result$Gene, data$misc$gene_annotations$input), "name"
        ]

        res$set_body(result)
    }
)

backend <- BackendRserve$new()
backend$start(app, http_port = 4000)
