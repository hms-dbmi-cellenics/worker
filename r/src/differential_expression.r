
# bh.adjust function from [1]
#' @description BH P-value adjustment with a log option. 
#' Helper functions for pagoda2 DE pipeline. 
bh.adjust <- function (x, log = FALSE) {
    nai <- which(!is.na(x))
    ox <- x
    x <- x[nai]
    id <- order(x, decreasing = FALSE)
    if (log) {
        q <- x[id] + log(length(x)/seq_along(x))
    }
    else {
        q <- x[id] * length(x)/seq_along(x)
    }
    a <- rev(cummin(rev(q)))[order(id)]
    ox[nai] <- a
    return(ox)
}

# compute_Zscore function from [1]
#' @description DE pipeline replicate from pagoda2.  
#' @param counts Counts matrix with cells as rows and genes as columns
#' @param groups gruops of cells to compare. Use the names of the cells
#' @param z.threshold numeric Minimal absolute Z score (adjusted) to report (default=3)
#' @param verbose boolean Whether to give verbose output (default=FALSE)
#' @param upregulated.only boolean Whether to report only genes that are expressed significantly higher in each group (default=FALSE)
#' 
#' @return Dataframe with columns:
#'              zscore
#'              log2fc
#'              pct
#'              Gene
#'              abszscore
#'              qval
compute_Zscore <- function(counts, groups, z.threshold=3, verbose = F, upregulated.only=FALSE){
    
    # Seurat give cells as columns, hence we need to transpose
    cm <- Matrix::t(counts) 
    
    if (!all(rownames(cm) %in% names(groups))) { 
        warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")
    }
    # determine a subset of cells that's in the groups and groups[cell]!=NA
    valid.cells <- rownames(cm) %in% names(groups)[!is.na(groups)]
    
    if (!all(valid.cells)) {
        # take a subset of the count matrix
        cm <- cm[valid.cells, ]
    }
    # reorder groups
    groups <- as.factor(groups[match(rownames(cm),names(groups))])
    
    groups <- as.factor(groups)
    if (verbose) {
        message("running differential expression with ",length(levels(groups))," clusters ... ")
    }
    
    # run wilcoxon test comparing each group with the rest
    lower.lpv.limit <- -100
    # calculate rank per-column (per-gene) average rank matrix
    n <- base::diff(cm@p)  ## number of non-zeros per column
    lst <- split(cm@x, rep.int(1:ncol(cm), n))  ## columns to list
    r_v2 <- unlist(lapply(lst, rank))  ## column-wise ranking and result collapsing
    xr <- cm  ## copy sparse matrix
    xr@x <- r_v2  ## replace non-zero elements with rank
    # calculate rank sums per group
    grs <- colSumByFac(xr,as.integer(groups))[-1,,drop=FALSE]
    # calculate number of non-zero entries per group
    xr@x <- numeric(length(xr@x))+1
    gnzz <- colSumByFac(xr,as.integer(groups))[-1,,drop=FALSE]
    #group.size <- as.numeric(tapply(groups,groups,length));
    group.size <- as.numeric(tapply(groups,groups,length))[1:nrow(gnzz)]
    group.size[is.na(group.size)]<-0 # trailing empty levels are cut off by colSumByFac
    # add contribution of zero entries to the grs
    gnz <- (group.size-gnzz)
    # rank of a 0 entry for each gene
    zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
    ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2
    # standardize
    n1n2 <- group.size*(nrow(cm)-group.size)
    # usigma <- sqrt(n1n2*(nrow(cm)+1)/12) # without tie correction
    # correcting for 0 ties, of which there are plenty
    usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
    usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
    x <- t((ustat - n1n2/2)/usigma) # standardized U value- z score
    
    # correct for multiple hypothesis
    if (verbose) {
        message("adjusting p-values ... ")
    }
    x <- matrix(qnorm(bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),ncol=ncol(x))*sign(x)
    rownames(x) <- colnames(cm)
    colnames(x) <- levels(groups)[1:ncol(x)]
    if (verbose) {
        message("done.\n")
    }
    
    # add fold change information
    log.gene.av <- log2(Matrix::colMeans(cm))
    group.gene.av <- colSumByFac(cm,as.integer(groups))[-1,,drop=FALSE] / (group.size+1)
    log2.fold.change <- log2(t(group.gene.av)) - log.gene.av
    # fraction of cells expressing
    f.expressing <- t(gnzz / group.size)
    max.group <- max.col(log2.fold.change)
    
    ds <- lapply(1:ncol(x),function(i) {
        z <- x[,i]
        vi <- which((if (upregulated.only) z else abs(z)) >= z.threshold)
        r <- data.frame(zscore=z[vi],log2fc=log2.fold.change[vi,i],pct=f.expressing[vi,i], Gene=rownames(x)[vi])
        rownames(r) <- r$Gene
        r$abszscore <- abs(r$zscore)
        
        r["qval"] <- lapply(
            r["abszscore"], function(x) 2 * pnorm(-x)
        )
        
        r["qval"] <- apply(
            r["qval"], 1, function(x) format(x, scientific = TRUE)
        )
        
        r <- r[order(r$zscore,decreasing=TRUE), ]
        r
    })
    names(ds) <- colnames(x)
    
    return(ds)
}

# runDE function
#' @description DE pipeline replicate from pagoda2.  
#' @param counts Counts matrix with cells as rows and genes as columns
#' @param groups gruops of cells to compare. Use the names of the cells
#' @param z.threshold numeric Minimal absolute Z score (adjusted) to report (default=3)
#' @param verbose boolean Whether to give verbose output (default=FALSE)
#' @param upregulated.only boolean Whether to report only genes that are expressed significantly higher in each group (default=FALSE)
#' 
#' @return Dataframe with columns:
#'              zscore
#'              log2fc
#'              pct
#'              Gene
#'              abszscore
#'              qval
runDE <- function(req){
    
    # set up a factor with the appropriate cells
    factor_DE <- rep(
        c("base", "background"),
        c(
            length(req$body$baseCells),
            length(req$body$backgroundCells)
        )
    ) 

    names(factor_DE) <- c( 
        req$body$baseCells, req$body$backgroundCells
    )

    # check mapped cells
    message("Mapped cells ", round(sum(names(factor_DE)%in%colnames(data))/ncol(data)*100, 4), " %.")

    # add seurat object the new groups to compare    
    data@meta.data$custom <- NA
    data@meta.data[names(factor_DE), "custom"] <- factor_DE
    
    # compute differential expression
    result <- compute_Zscore(data@assays$RNA@data, groups = factor_DE)$base

    # replace name with Gene names
    result$gene_names <- data@misc$gene_annotations[
        match(rownames(result), data@misc$gene_annotations$input), "name"
    ]
    
    # Alternative with seurat DE pipeline
    # result <- FindMarkers(data, group.by = "custom", ident.1 = "base", ident.2 = "background")
    # # replace name with Gene names
    # result$gene_names <- data@misc$gene_annotations[
    #    match(rownames(result), data@misc$gene_annotations$input), "name"
    #]
    # result$Gene <- rownames(result)
    # result$p_val <- NULL
    # colnames(result) <- c("log2fc", "pct.1", "pct.2", "p_val_adj", "gene_names", "Gene")

    return(result)
}


# [1]   Replicated from package pagoda2
#       Nikolas Barkas, Viktor Petukhov, Peter Kharchenko and Evan
#       Biederstedt (2020). pagoda2: Single Cell Analysis and Differential
#       Expression. R package version 1.0.0.
#       https://github.com/kharchenkolab/pagoda2