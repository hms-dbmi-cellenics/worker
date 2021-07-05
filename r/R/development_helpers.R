generate_request <- function(task) {
    res = switch(task,
           "expression" = list(body=list(genes=list("gzma")))
           )
}

load_experiment_data <- function(exp_id="test"){
    return(readRDS(paste0("../data/",exp_id,"/r.rds",sep="")))
}
