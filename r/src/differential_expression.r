
# runDE function
# Differential expression analysis can be requested. The request post include the cells ID that we need
# to compare in a cellSet named "base" and "background",
# For now, we are going to support Wilcoxon Rank Sum test, but we can add later other tests like "t", "MAST"
# or "DESeq2". We are going to return all the genes that have a absolute value avg_log2FC greater than 0.25 (this
# default value becomes from the Seurat default value of function FindMarkers)
# Out Gene, avg_log2FC,  p_val_adj,  pct_1, pct_2
#' @description DE pipeline replicate
#' 
#' @return Dataframe with columns:
#'              Gene
#'              avg_log2FC
#'              p_val_adj
#'              pct_1
#'              pct_2
#'              zscore          <-- Need to be returned until the UI is changed
#'              log2fc          <-- Need to be returned until the UI is changed
#'              pct             <-- Need to be returned until the UI is changed
#'              abszscore       <-- Need to be returned until the UI is changed
#'              qval            <-- Need to be returned until the UI is changed
runDE <- function(req){
    
    cells_id <- data@meta.data$cells_id

    # Remove filtered cells
    baseCells <- req$body$baseCells[req$body$baseCells %in% cells_id]
    backgroundCells <- req$body$backgroundCells[req$body$backgroundCells %in% cells_id]

    # set up a factor with the appropriate cells
    factor_DE <- rep(
        c("base", "background"),
        c(
            length(baseCells),
            length(backgroundCells)
        )
    ) 


    baseCells_barcode <- rownames(data@meta.data)[match(baseCells, cells_id)]
    backgroundCells_barcode <- rownames(data@meta.data)[match(backgroundCells, cells_id)]

    names(factor_DE) <- c( 
        baseCells_barcode, backgroundCells_barcode
    )

    # check mapped cells
    message("Mapped cells ", round(sum(names(factor_DE)%in%colnames(data))/ncol(data)*100, 4), " %.")

    # add seurat object the new groups to compare    
    data@meta.data$custom <- NA
    data@meta.data[names(factor_DE), "custom"] <- factor_DE

    # Compute differential expression
    result <- FindMarkers(data, group.by = "custom", ident.1 = "base", ident.2 = "background")

    # Replace name with Gene names
    result$gene_names <- data@misc$gene_annotations[
        match(rownames(result), data@misc$gene_annotations$input), "name"
    ]
    result$Gene <- rownames(result)

    # As a first view, order by p_val_adj, to have the most significant at first.
    result <- result[order(result$p_val_adj, decreasing = F), ]

    # Change "." in pct.1 by _
    colnames(result) <- gsub("[.]", "_", colnames(result))

    # Check if the gene_symbol does not appear in annotation. In that case the NA value will be changed to ENSEMBL ID
    result$gene_names[is.na(result$gene_names)] <- result$Gene[is.na(result$gene_names)]


    ## Old DE results from pagoda2
    #result$zscore <- result$pct_1
    #result$pct <- result$pct_1
    #result$abszscore <- result$pct_2
    #result$log2fc <- result$avg_log2FC
    #result$qval <- result$p_val_adj

    return(result)
}
