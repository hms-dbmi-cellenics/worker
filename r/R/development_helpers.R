#' Generates request for worker task
#'
#' @param task
#'
#' @return
#' @export
#'
#' @examples
generate_request <- function(task) {
  res <- switch(task,
    "expression" = list(body = list(genes = list("gzma", "Lyz2"))),
    "markerHeatmap" = list(body = list(nGenes = 5, type = "louvain", config = list(resolution = 0.5))),
    "DE" = list(body = list(baseCells = 0:500, backgroundCells = 501:1499))
  )
}

#' Loads Experiment data
#'
#' @param exp_id
#'
#' @return
#' @export
#'
#' @examples
load_experiment_data <- function(exp_id = "test") {
  return(readRDS(paste0("../data/", exp_id, "/r.rds", sep = "")))
}
