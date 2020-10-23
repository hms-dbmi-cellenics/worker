library(RestRserve)
library(pagoda2)
library(conos)
library(Matrix)
require(data.table)
library(DESeq2)
library(RJSONIO)

load_data <- function() {
    experiment_id <- Sys.getenv("EXPERIMENT_ID")
    message(paste("Welcome to Biomage R worker, experiment id", experiment_id))

    loaded <- F
    data <- NULL

    while (!loaded) {
        data <- tryCatch(
            {
                f <- readRDS(
                    paste(
                        "/data", experiment_id, "r.rds",
                        sep = "/"
                    )
                )

                loaded <- T
                length <- dim(f$counts)

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

            result <- result[c("Z", "M", "Gene")]
            colnames(result) <- c("zscore", "log2fc", "gene_names")

            # compute absolute Z score from Z score
            result["abszscore"] <- abs(result["zscore"])

            # compute q-value from absolute z score
            result["qval"] <- lapply(
                result["abszscore"], function(x) 2 * pnorm(-x)
            )

            result["qval"] <- apply(
                result["qval"], 1, function(x) format(x, scientific = TRUE)
            )

            message("yolo")

            res$set_body(result)
        }
    )

    return(app)
}

# apply(hist_data, 2, function(x) pnorm(x, mean=mean(x), sd=sd(x)))

data <- load_data()
backend <- BackendRserve$new()

backend$start(create_app(data), http_port = 4000)