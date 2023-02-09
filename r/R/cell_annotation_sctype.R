library("openxlsx")
library("HGNChelper")

format_sctype_cell_sets_object <-
  function(data, species, tissue, color_pool) {
    cell_class_key <- paste0("ScType-", tissue, "-", species)

    cell_class <-
      list(
        key = cell_class_key,
        name = cell_class_key,
        rootNode = TRUE,
        type = "cellSets",
        children = list()
      )

    for (i in 1:length(unique(data@meta.data$customclassif))) {
      cell_set_key <- unique(data@meta.data$customclassif)[[i]]

      new_cell_set <- list(
        key = cell_set_key,
        name = cell_set_key,
        rootNode = FALSE,
        type = "cellSets",
        color = sample(data@misc$color_pool, 1), 
        cellIds = data@meta.data[data@meta.data$customclassif == data@meta.data$customclassif[[i]], "cells_id"]
      )
      color_pool <- color_pool[-1]
      cell_class$children <- append(cell_class$children, list(new_cell_set))
    }

    return(cell_class)
  }

ScTypeAnnotate <- function(req, data) {
  cell_sets <- req$body$cellSets
  species <- req$body$species
  tissue <- req$body$tissue

  if ("integrated" %in% names(data@assays)) {
    active_assay <- "integrated"
  } else if ("SCT" %in% names(data@assays)) {
    active_assay <- "SCT"
  } else {
    active_assay <- "RNA"
  }

  scale_data <- get_formatted_data(data, active_assay)
  parsed_cellsets <- parse_cellsets(cell_sets)
  data <- add_clusters(data, parsed_cellsets)

  data[[active_assay]]@scale.data <- scale_data

  data <- run_sctype(data, active_assay, tissue, species)

  message("formatting cellsets")
  formatted_cell_class <- format_sctype_cell_sets_object(data, species, tissue, data@misc$color_pool)

  message("updating through api")
  updateCellSetsThroughApi(
    formatted_cell_class,
    req$body$apiUrl,
    data@misc$experimentId,
    formatted_cell_class$key,
    req$body$authJwt
  )

  return(formatted_cell_class)
}


get_formatted_data <- function(scdata, active_assay) {
  scale_data <- data.table::as.data.table(scdata[[active_assay]]@scale.data, keep.rownames = TRUE)

  # convert ENS ID to gene symbol
  scale_data <- add_gene_symbols(scale_data, scdata)

  # take gene with highest mean expression for duplicated gene symbols
  scale_data <- collapse_genes(scale_data)

  # convert to matrix and set gene symbols as rownames of the count matrix
  scale_data <- format_matrix(scale_data)

  return(scale_data)
}


add_gene_symbols <- function(scale_data, scdata) {
  colnames(scale_data)[1] <- "input"
  annot <- data.table::setDT(scdata@misc$gene_annotations)
  annot <- annot[, .(input, original_name)]

  feature_types <- get_feature_types(annot)
  if (feature_types == IDS_IDS) {
    stop(
      generateErrorMessage(
        error_codes$NO_GENE_SYMBOLS,
        "Features file doesn't contain gene symbols."
      )
    )
  }

  if (feature_types == SYM_IDS) {
    annot[, c(1, 2)] <- annot[, c(2, 1)]
  }

  scale_data <- annot[scale_data, on = .(input)]

  return(scale_data)
}

collapse_genes <- function(scale_data) {
  scale_data[, mean_expression := rowMeans(.SD), .SDcols = !c("input", "original_name")]
  scale_data <- scale_data[!duplicated(scale_data, by = "original_name") | mean_expression == max(mean_expression), ]
  scale_data[, mean_expression := NULL]

  return(scale_data)
}


format_matrix <- function(scale_data) {
  gene_symbols <- scale_data[, original_name]
  scale_data <- apply(as.matrix.noquote(scale_data[, -c(1, 2)]), 2, as.numeric)
  rownames(scale_data) <- gene_symbols

  return(scale_data)
}


run_sctype <- function(data, active_assay, tissue, species) {
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  db <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  gs_list <- gene_sets_prepare(db, tissue)

  # get cell-type by cell matrix
  cell_type_scores <- sctype_score(
    scRNAseqData = data[[active_assay]]@scale.data, scaled = TRUE,
    gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
  )

  # merge by cluster
  clusters <- "seurat_clusters"
  metadata_clusters <- data@meta.data[[clusters]]
  cluster_scores <- do.call("rbind", lapply(unique(metadata_clusters), function(cl) {
    cell_type_scores_cl <- sort(rowSums(cell_type_scores[, as.integer(rownames(data@meta.data[metadata_clusters == cl, ]))]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(cell_type_scores_cl), scores = cell_type_scores_cl, ncells = sum(metadata_clusters == cl)), 10)
  }))

  sctype_scores <- cluster_scores %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)


  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"

  data@meta.data$customclassif <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    data@meta.data$customclassif[metadata_clusters == j] <- as.character(cl_type$type[1])
  }

  return(data)
}
