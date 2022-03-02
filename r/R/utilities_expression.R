#' Returns expression values for selected genes
#'
#' @param genes - Must have names and input(ensmbl ids)
#' @param data
#'
#' @return
#' @export
#'
#' @examples
getExpressionValues <- function(genes, data) {
  quantile_threshold <- 0.95
  #
  # Get the expression values for those genes in the corresponding matrix.
  geneExpression <- data@assays$RNA@data[unique(genes$input), , drop = FALSE]
  geneExpression <- as.data.frame(t(as.matrix(geneExpression)))
  geneExpression$cells_id <- data@meta.data$cells_id
  geneExpression <- geneExpression[order(geneExpression$cells_id), ]
  geneExpression <- geneExpression %>%
    tidyr::complete(cells_id = seq(0, max(data@meta.data$cells_id))) %>%
    dplyr::select(-cells_id) %>%
    as.data.frame()
  # worried about duplicate gene row.names in @data
  symbol_idx <- match(colnames(geneExpression), genes$input)
  colnames(geneExpression) <- genes$name[symbol_idx]
  adjGeneExpression <- as.data.frame(apply(geneExpression, 2, FUN = function(x) {
    lim <- as.numeric(quantile(x, quantile_threshold, na.rm = TRUE))
    i <- 0.01
    while (lim == 0 & i + quantile_threshold <= 1) {
      lim <- as.numeric(quantile(x, quantile_threshold + i, na.rm = TRUE))
      i <- i + 0.01
    }
    return(pmin(x, lim))
  }))
  return(list(rawExpression = geneExpression, truncatedExpression = adjGeneExpression))
}

#' formatExpression
#'
#' @param data
#'
#' @return Formats expression values for the UI.
#' @export
#'
#' @examples
formatExpression <- function(data) {
  return(purrr::map2(data$rawExpression, data$truncatedExpression, formatExpressionAux))
}

formatExpressionAux <- function(raw, trunc) {
  return(list(
    rawExpression = list(
      mean = mean(raw, na.rm = TRUE),
      stdev = sd(raw, na.rm = TRUE),
      expression = raw
    ),
    truncatedExpression = list(
      min = min(trunc, na.rm = TRUE),
      max = max(trunc, na.rm = TRUE),
      expression = trunc
    )
  ))
}
