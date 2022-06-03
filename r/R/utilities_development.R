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
