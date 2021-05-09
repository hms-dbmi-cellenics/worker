
# Function to retrieve all the doublet score for the current experiment.
# The doublet scores were computing in the data-ingest script. To compute then, we
# use the package Scrublet [1]
getDoubletScore <- function(req) {
  # Check if the experiment has doublet_scores stored in rds file
  if(!"doublet_scores"%in%colnames(data@meta.data)){
    stop("Doublet scores are not computed for this experiment.")
  }
  
  # Subset the doublet_scores ordering by cells_id
  result <- data@meta.data[order(data$cells_id, decreasing = F), "doublet_scores"]
  result<- as.data.frame(result)
  result$cells_id <- data@meta.data$cells_id
  result <- result %>% complete(cells_id = seq(0,max(data@meta.data$cells_id))) %>% select(-cells_id)
  result <- t(unname(result))
  result <- purrr::map(result, function(x) { if(is.na(x)) { return(NULL) } else { return(x) } } )
  
  # Be aware of possible na values
  if(any(is.null(result)))
    warning("There are missing values in the doublet_scores results")
  return(result)
  
}

# Function to retrieve all the mt-content score for the current experiment.
# The MT-content was computing in the data-ingest script. To compute then, we
# use the function PercentageFeatureSet fom Seurat package. We have been able to identify
# the MT-genes only in MMusculus and Homo Sapiens by grepping "MT"
getMitochondrialContent <- function(req) {
    
    # Check if the experiment has percent.mt stored in rds file
    if(!"percent.mt"%in%colnames(data@meta.data)){
        stop("MT content is not computed for this experiment")
    }

    # Subset the percent.mt ordering by cells_id (DESC)
    result <- data@meta.data[order(data$cells_id, decreasing = F), "percent.mt"]


    # Be aware of possible na values
    if(any(is.na(result)))
        warning("There are missing values in the doublet_scores results")
    
    return(result)

}


# [1] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell 
# Transcriptomic Data. Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3. 
# PMID: 30954476; PMCID: PMC6625319.

