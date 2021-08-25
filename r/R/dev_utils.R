# see what scaled heatmap should look like
#
# scale.data: output of getHeatmapExpression
# clusters: data$seurat_clusters
#
plotScaledHeatmap <- function(scale.data, clusters) {
    cell.order <- unique(names(sort(clusters)))
    gene.order <- rev(row.names(scale.data))

    scale.data <- scale.data[complete.cases(scale.data), ]

    scale.data <- scale.data %>%
        as_tibble(rownames = 'feature') %>%
        tidyr::pivot_longer(-feature, names_to = 'cell', values_to = 'expression') %>%
        mutate('identity' = clusters[cell],
               'cell' = factor(cell, levels = cell.order),
               'feature' = factor(feature, levels = gene.order))

    limits <- c(min(scale.data$expression), max(scale.data$expression))
    colors <- Seurat::PurpleAndYellow()

    ggplot(scale.data) +
        geom_raster(mapping = aes_string(x = "cell", y = "feature", fill = "expression")) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_fill_gradientn(limits = limits,
                             colors = colors,
                             na.value = "white") +
        labs(x = NULL, y = NULL, fill = "expression") +
        Seurat::WhiteBackground() +
        Seurat::NoAxes(keep.text = FALSE)
}
