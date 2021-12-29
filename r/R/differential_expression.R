#' runDE function
#'
#' Differential expression analysis can be requested. The request post include the cells ID that we need
#' to compare in a cellSet named "base" and "background",
#' For now, we are going to support Wilcoxon Rank Sum test, but we can add later other tests like "t", "MAST"
#' or "DESeq2". We are going to return all the genes that have a absolute value avg_log2FC greater than 0.25 (this
#' default value becomes from the Seurat default value of function FindMarkers)
#' Out Gene, avg_log2FC,  p_val_adj,  pct_1, pct_2
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
#' @export
#'
runDE <- function(req, data) {

  # add comparison group to 'custom' slot
  data <- addComparisonGroup(req, data)

  # Compute differential expression
  result <- presto::wilcoxauc(data, assay = "data", seurat_assay = "RNA", group_by = "custom")
  result <- result[result$group == "base", ]

  rownames(result) <- result$feature
  result <- result[, c("pval", "logFC", "pct_in", "pct_out", "padj", "auc")]
  colnames(result) <- list("p_val", "logFC", "pct_1", "pct_2", "p_val_adj", "auc")

  message("checking FindMarkers results:  ", str(result))
  # Replace name with Gene names
  result$gene_names <- data@misc$gene_annotations[
    match(rownames(result), data@misc$gene_annotations$input), "name"
  ]
  result$Gene <- rownames(result)

  # As a first view, order by p_val_adj, to have the most significant at first.
  result <- result[order(result$p_val_adj, decreasing = FALSE), ]

  # Change "." in pct.1 by _
  colnames(result) <- gsub("[.]", "_", colnames(result))

  # Check if the gene_symbol does not appear in annotation. In that case the NA value will be changed to ENSEMBL ID
  result$gene_names[is.na(result$gene_names)] <- result$Gene[is.na(result$gene_names)]


  if (!("pagination" %in% names(req$body))) {
    result <- list(gene_results = purrr::transpose(result),full_count=nrow(result))
    message("Pagination not enabled, returning results: ", str(result))
    return(result)
  }

  message("Paginating results:  ", str(result))
  pagination <- req$body$pagination
  pathwayAnalysis <- req$body$pathwayAnalysis
  order_by <- pagination$orderBy
  order_decreasing <- pagination$orderDirection == "DESC"
  offset <- pagination$offset
  limit <- pagination$limit
  filters <- pagination$filters

  result <- applyFilters(result, filters)

  if("pathwayAnalysis" %in% names(req$body)){
    if(pathwayAnalysis$returnOnlyGenes){
      n_genes = pathwayAnalysis$nGenesToReturn
      result <- list(gene_results=result$Gene[1:n_genes],full_count=n_genes)
    }
  }else{
    result <- handlePagination(result, offset, limit, order_by, order_decreasing)
    result$gene_results <- purrr::transpose(result$gene_results)
  }

  return(result)
}
