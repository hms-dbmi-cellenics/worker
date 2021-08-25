
# modeled after Seurat::DoHeatmap
getHeatmapExpression <- function(enids, data, disp.min = -2.5, disp.max = 2.5) {
    data <- data[enids, ]
    data <- Seurat::ScaleData(data, features = row.names(data))

    scale.data <- data[['RNA']]@scale.data %>%
        Seurat::MinMax(disp.min, disp.max)

    return(scale.data)
}

# adds NA rows for filtered cells
completeCellIds <- function(expr, cells_id) {
    expr <- as.data.frame(t(expr))
    row.names(expr) <- cells_id
    idx <- as.character(0:max(cells_id))
    expr <- expr[idx,, drop = FALSE]
    row.names(expr) <- idx
    return(expr)
}

# caps rows of expression matrix (genes) to quantile threshold
truncateExpression <- function(expr, quantile_threshold = 0.95) {

    adj <- apply(expr, 1, FUN = function(x) {
        lim <- quantile(x, quantile_threshold, na.rm = TRUE)
        i <- 0.01
        while (lim == 0 & i + quantile_threshold <= 1) {
            lim <- quantile(x, quantile_threshold + i, na.rm = TRUE)
            i <- i + 0.01
        }
        return(pmin(x, lim))
    })
    return(t(adj))
}
