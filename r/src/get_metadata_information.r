getDoubletScore <- function(req) {
    
    cells <- req$body$cells
    
    if(any(!cells%in%rownames(data@meta.data))){
        stop("There are some requested cells that are not in the data")
    }

    if(!"doublet_scores"%in%colnames(data@meta.data)){
        stop("Doublet scores are not computed for this experiment")
    }

    result <- subset(data@meta.data, rownames(data@meta.data)%in%cells, "doublet_scores")

    if(any(is.na(result[, "doublet_scores"])))
        warning("There are missing values in the doublet_scores results")

    return(result)

}


getMitochondrialContent <- function(req) {
    
    cells <- req$body$cells
    
    if(any(!cells%in%rownames(data@meta.data))){
        stop("There are some requested cells that are not in the data")
    }

    if(!"percent.mt"%in%colnames(data@meta.data)){
        stop("MT content is not computed for this experiment")
    }

    result <- subset(data@meta.data, rownames(data@meta.data)%in%cells, "percent.mt")
    
    if(any(is.na(result[, "percent.mt"])))
        warning("There are missing values in the doublet_scores results")

    
    return(result)

}
