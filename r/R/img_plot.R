runRidgePlot <- function(req,data,config){
  data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = "pca")
  data <- Seurat::FindClusters(data, resolution = 0.8, verbose = FALSE, algorithm = "louvain")

  df <- data@misc$gene_annotations

  genesSubset <- subset(df, toupper(df$name) %in% toupper(req$genes))

  if (!nrow(genesSubset)) {
    stop(
      generateErrorMessage(
        error_codes$GENE_NOT_FOUND,
        paste("Gene(s):", paste(req$genes, collapse = ", "), "not found!")
      )
    )
  }

  genesSubset <- genesSubset[,"input"]
  Seurat::RidgePlot(data,genesSubset)
  ggplot2::ggsave("./plot.png")

  raw_img <- loder::readPng("./plot.png")

  return(RJSONIO::toJSON(list(data=raw_img)))
}

runVlnPlot <- function(req,data,config){
  data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = "pca")
  data <- Seurat::FindClusters(data, resolution = 0.8, verbose = FALSE, algorithm = "louvain")

  df <- data@misc$gene_annotations

  genesSubset <- subset(df, toupper(df$name) %in% toupper(req$genes))

  if (!nrow(genesSubset)) {
    stop(
      generateErrorMessage(
        error_codes$GENE_NOT_FOUND,
        paste("Gene(s):", paste(req$genes, collapse = ", "), "not found!")
      )
    )
  }

  genesSubset <- genesSubset[,"input"]
  Seurat::VlnPlot(data,genesSubset)
  ggplot2::ggsave("./plot.png")

  raw_img <- loder::readPng("./plot.png")

  return(RJSONIO::toJSON(list(data=raw_img)))
}

runDotPlot <- function(req,data,config){
  data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = "pca")
  data <- Seurat::FindClusters(data, resolution = 0.8, verbose = FALSE, algorithm = "louvain")

  aucMin = 0.3
  pctInMin = 20
  pctOutMax = 70
  nFeatures <- 2

  all_markers <- presto::wilcoxauc(data, group_by = "seurat_clusters", assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)

  # may not return nFeatures markers per cluster if values are too stringent
  filtered_markers <- all_markers %>%
    dplyr::filter(logFC > 0 &
                    auc >= aucMin &
                    pct_in >= pctInMin &
                    pct_out <= pctOutMax) %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(pval))

  top_markers <- filtered_markers %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::arrange(group)


  Seurat::DotPlot(data,features=top_markers$feature,)
  ggplot2::ggsave("./plot.png")

  raw_img <- loder::readPng("./plot.png")

  return(RJSONIO::toJSON(list(data=raw_img)))
}

runMarkerHeat <- function(req,data,config){
  data <- Seurat::FindNeighbors(data, annoy.metric = "cosine", verbose = FALSE, reduction = "pca")
  data <- Seurat::FindClusters(data, resolution = 0.8, verbose = FALSE, algorithm = "louvain")

  aucMin = 0.3
  pctInMin = 20
  pctOutMax = 70
  nFeatures <- 2

  all_markers <- presto::wilcoxauc(data, group_by = "seurat_clusters", assay = "data", seurat_assay = "RNA")
  all_markers$group <- as.numeric(all_markers$group)

  # may not return nFeatures markers per cluster if values are too stringent
  filtered_markers <- all_markers %>%
    dplyr::filter(logFC > 0 &
                    auc >= aucMin &
                    pct_in >= pctInMin &
                    pct_out <= pctOutMax) %>%
    dplyr::group_by(feature) %>%
    dplyr::slice(which.min(pval))

  top_markers <- filtered_markers %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = nFeatures) %>%
    dplyr::arrange(group)


  Seurat::DoHeatmap(subset(data, downsample = 100), features = top_markers$feature, size = 3)
  ggplot2::ggsave("./plot.png")

  raw_img <- loder::readPng("./plot.png")

  return(RJSONIO::toJSON(list(data=raw_img)))
}


runFeaturePlot <- function(req,data,config){
  df <- data@misc$gene_annotations
  cell_sets <- req$cell_sets

  for (cellset in cell_sets) {
    if (cellset$key == "louvain") { 
      clusters <- cellset$children
    }
  }
  cell_ids <- c()
  name <- c()

  for(cluster in clusters) {
    cell_ids <- c(cell_ids, cluster$cellIds)
    name <- c(name, rep(cluster$name,length(cluster$cellIds)))
  }


  cluster_table <- data.frame(cell_ids, name)
  print(str(data))
  print(tail(cluster_table))
  genesSubset <- subset(df, toupper(df$name) %in% toupper(req$genes))

  if (!nrow(genesSubset)) {
    stop(
      generateErrorMessage(
        error_codes$GENE_NOT_FOUND,
        paste("Gene(s):", paste(req$genes, collapse = ", "), "not found!")
      )
    )
  }

  genesSubset <- genesSubset[,"input"]
  Seurat::FeaturePlot(data,genesSubset)
  ggplot2::ggsave("./plot.png")


  raw_img <- loder::readPng("./plot.png")
  return(RJSONIO::toJSON(list(data=raw_img)))

}

createImgPlot <- function(req, data) {
  config = list(
    credentials = list(
      creds = list(
        access_key_id = "mock-access-key",
        secret_access_key = "mock-secure-acces-key"
      ),
      profile = "string"
    ),
    endpoint = "http://host.docker.internal:4566",
    region = "eu-west-1")
  req <- req$body
  imgPlotTypes <- list("ridgePlot"=runRidgePlot,"featurePlot"=runFeaturePlot,"dotPlot"=runDotPlot,"vlnPlot"=runVlnPlot,"markerHeatmap"=runMarkerHeat)

  res<-imgPlotTypes[[req$plotSubType]](req,data,config)
  return(res)
}

# put_object_in_s3 <- function(config, bucket, object, key,data) {
#   message(sprintf("Putting %s in %s", key, bucket))
#   message(config)
#   s3 <- paws::s3(config = config)
#   s3$put_object(
#     Bucket = bucket,
#     Key = key,
#     Body = object
#   )
#   tag <- list(TagSet=list(
#     list("Key"="experimentId","Value"=data@misc$experimentId),
#     list("Key"="requestType","Value"="ImgPlot"))
#   )
#   print(tag)
#   s3$put_object_tagging(bucket, key,Tagging=tag)
# }


