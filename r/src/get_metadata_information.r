getDoubletScore <- function(req) {
     
    # Check cells
    cells <- check_valid_cells(req$body$cells)

    # Check if the experiment has doublet_scores stored in rds file
    if(!"doublet_scores"%in%colnames(data@meta.data)){
        stop("Doublet scores are not computed for this experiment.")
    }

    # Subset the doublet_scores for the requested cells
    result <- subset(data@meta.data, data$cells_id%in%cells, c("cells_id", "doublet_scores"))

    # Be aware of possible na values
    if(any(is.na(result[, "doublet_scores"])))
        warning("There are missing values in the doublet_scores results")

    return(result)

}


getMitochondrialContent <- function(req) {
    
    # Check cells
    cells <- check_valid_cells(req$body$cells)

    # Check if the experiment has percent.mt stored in rds file
    if(!"percent.mt"%in%colnames(data@meta.data)){
        stop("MT content is not computed for this experiment")
    }

    # Subset the percent.mt for the requested cells
    result <- subset(data@meta.data, data$cells_id%in%cells, c("cells_id", "percent.mt"))
    
    # Be aware of possible na values
    if(any(is.na(result[, "percent.mt"])))
        warning("There are missing values in the doublet_scores results")
    
    return(result)

}


# Internal functio to check wheter the cells belongs to the experiment and if the experiment 
# can retrieve information from them.
check_valid_cells <- function(cells){
    
    # Convert to numeric if required
    if(!is.numeric(cells))
        cells <- as.numeric(cells)

    # Check any possible NA value
    if(any(is.na(cells))){
        cells <- na.omit(cells)
        warning("There is NA values in the request.")
    }

    # Check if the experiment has cells_id stored in rds file
    if(!"cells_id"%in%colnames(data@meta.data)){
        stop("Cells id are not computed for this experiment.")
    }

    # Check cells that are not in the experiment
    if(any(!cells%in%data$cells_id)){
        stop("There are some requested cells that are not in the data: ", cells[!cells%in%data$cells_id])
    }

    return(cells)
}
