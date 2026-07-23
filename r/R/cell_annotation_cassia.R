#' Annotate cells with CASSIA
#'
#' This function uses [CASSIA](https://github.com/ElliotXie/CASSIA),
#' an LLM based collaborative agent system, to annotate cells in a
#' Seurat object.
#'
#' As a preprocessing step, it parses the cellsets object, adds the cluster
#' information as metadata columns and sets the Louvain clusters as the active
#' identity. It then computes the top marker genes per cluster (subsampling
#' cells for speed) and passes them to the CASSIA pipeline. The cluster-level
#' cell type predictions are added back to the metadata slot and used to build
#' a new cell class, which is finally patched into the cellsets object via
#' the API.
#'
#' It is modelled after [ScTypeAnnotate()], differing only in the annotation
#' engine (CASSIA instead of ScType) and in accepting free text species/tissue.
#'
#' @param req list with items:
#'  \code{req$body$cellSets} cellsets object
#'  \code{req$body$species} free text species (e.g. "Human", "Mouse")
#'  \code{req$body$tissue} free text tissue type (e.g. "Large Intestine")
#'
#' @param data Seurat object
#'
#' @export
#'
CASSIAAnnotate <- function(req, data) {
  cell_sets <- req$body$cellSets
  species <- req$body$species
  tissue <- req$body$tissue
  # Optional free-text study context (e.g. "3 colorectal tumor + 2 normal
  # adjacent tissue samples"), passed to the CASSIA agents to inform annotation.
  additional_info <- req$body$additionalInfo
  experiment_id <- req$body$experimentId
  task_name <- req$body$name

  # Bedrock model used for all CASSIA agents (annotation, scoring, boost, merge)
  model <- "qwen.qwen3-next-80b-a3b"

  # Phase progress is written to a shared file that the python worker polls to
  # emit UI heartbeats (see helpers/task_progress.py). The prep phases are fast
  # and are reported message-only (no percentage); only the CASSIA batch itself
  # reports a percentage, driven per-cluster by enable_cassia_batch_progress.
  progress_file <- worker_progress_file(experiment_id, task_name)
  on.exit(clear_worker_progress(experiment_id, task_name), add = TRUE)

  # Authenticate to Bedrock with the worker's IAM role (no static API key) and
  # get the OpenAI-compatible Bedrock endpoint to use as the CASSIA provider.
  provider <- set_cassia_bedrock_auth()

  # add cluster info to the metadata and set louvain clusters as active ident
  write_worker_progress(experiment_id, task_name, "Preparing clusters")
  children_cell_sets <- sapply(cell_sets, `[[`, "children")
  parsed_cellsets <- parse_cellsets(children_cell_sets)
  data <- add_clusters(data, parsed_cellsets, cell_sets)

  # Number the clusters 1..n and annotate on those numbers, so CASSIA predicts
  # cell types from marker genes alone rather than echoing the user's existing
  # cluster names (e.g. "Epithelial").
  cluster_levels <- unique(
    data@meta.data$seurat_clusters[!is.na(data@meta.data$seurat_clusters)]
  )
  data@meta.data$cassia_cluster <- match(
    data@meta.data$seurat_clusters, cluster_levels
  )

  # compute per-cluster marker genes and run the CASSIA annotation pipeline
  write_worker_progress(experiment_id, task_name, "Computing marker genes")
  markers <- get_cassia_markers(data, cluster_col = "cassia_cluster")

  write_worker_progress(
    experiment_id, task_name,
    "Running CASSIA annotation (this can take several minutes)"
  )
  data <- run_cassia(
    data, markers, tissue, species, model, provider,
    additional_info = additional_info,
    cluster_col = "cassia_cluster",
    progress_file = progress_file
  )

  # CASSIA produces up to three progressively coarser merged groupings; save
  # each as its own cell class (keyed CASSIA-<tissue>-<species>-1/-2/-3).
  write_worker_progress(
    experiment_id, task_name, "Saving annotated cell sets"
  )
  formatted_cell_classes <- list()
  for (grouping in 1:3) {
    annotation_column <- paste0("CASSIA_merged_grouping_", grouping)
    if (!annotation_column %in% colnames(data@meta.data)) next

    formatted_cell_class <- format_cassia_cell_sets(
      data, species, tissue, grouping
    )

    updateCellSetsThroughApi(
      formatted_cell_class,
      req$body$apiUrl,
      req$body$experimentId,
      formatted_cell_class$key,
      req$body$authJwt
    )

    formatted_cell_classes <- append(
      formatted_cell_classes, list(formatted_cell_class)
    )
  }
  message("CASSIAAnnotate: added ", length(formatted_cell_classes),
          " cell classes (merged groupings)")

  return(formatted_cell_classes)
}


#' Authenticate CASSIA to Amazon Bedrock using the worker's IAM role
#'
#' CASSIA calls Bedrock through its OpenAI-compatible endpoint
#' (`https://bedrock-runtime.<region>.amazonaws.com/openai/v1`), which
#' authenticates with a bearer token. Rather than a static API key, we mint a
#' short-term Bedrock bearer token from the worker's own IAM credentials (IRSA
#' role in-cluster, default AWS credential chain locally) using
#' `aws_bedrock_token_generator`. The token is generated client-side (no network
#' call) and is valid for up to 12h — far longer than a single annotation task.
#'
#' The token is written straight into Python's `os.environ` as
#' `CUSTOMIZED_API_KEY` (the var CASSIA's custom-HTTP provider reads). We set it
#' in Python directly because reticulate snapshots `os.environ` when it boots the
#' interpreter, so R's `Sys.setenv` would not propagate to CASSIA.
#'
#' @return the OpenAI-compatible Bedrock endpoint URL to pass as the CASSIA
#'  provider (`overall_provider`)
#' @export
#'
set_cassia_bedrock_auth <- function() {
  # CASSIA reuses the worker's existing reticulate virtualenv instead of its
  # default "cassia_env", so the long-lived R process binds a single Python.
  # Must be set before CASSIA's namespace loads (i.e. before any CASSIA:: call).
  options(CASSIA.env_name = "r-reticulate")

  # Treat an unset OR empty region env var as "use the default".
  region <- Sys.getenv("AWS_REGION", unset = "")
  if (!nzchar(region)) region <- Sys.getenv("AWS_DEFAULT_REGION", unset = "")
  if (!nzchar(region)) region <- "us-east-1"

  token <- tryCatch({
    gen <- reticulate::import("aws_bedrock_token_generator")
    gen$provide_token(region = region)
  }, error = function(e) {
    stop(generateErrorMessage(
      error_codes$CASSIA_ERROR,
      paste0(
        "Could not mint a Bedrock bearer token from the worker's IAM ",
        "credentials (region ", region, "). Check that the worker role has ",
        "Bedrock access and that AWS credentials are available: ",
        conditionMessage(e)
      )
    ))
  })

  os <- reticulate::import("os")
  os$environ$update(list(CUSTOMIZED_API_KEY = token))

  sprintf("https://bedrock-runtime.%s.amazonaws.com/openai/v1", region)
}


# Runtime hook (applied to CASSIA's in-memory BatchProgressTracker via
# reticulate) that mirrors per-cluster progress to CASSIA_PROGRESS_FILE as each
# cluster completes. We do NOT edit CASSIA's installed source (it may be
# installed from PyPI/CRAN, not our vendored copy) - the class method is wrapped
# at runtime. Idempotent, and every step is wrapped so a CASSIA internals change
# degrades gracefully to the coarser phase-level progress. The env var points at
# the generic per-task progress file (see worker_progress_file).
CASSIA_BATCH_PROGRESS_PATCH <- r"(
import os, sys, json, threading

def _cellenics_patch_batch_progress():
    def make_wrapped(orig):
        def wrapped(self, name):
            orig(self, name)
            path = os.environ.get("CASSIA_PROGRESS_FILE")
            if not path:
                return
            try:
                with self.lock:
                    payload = {
                        "step": self.completed,
                        "total": self.total,
                        "message": getattr(self, "title", "Annotating clusters"),
                    }
                tmp = path + "." + str(threading.get_ident()) + ".tmp"
                with open(tmp, "w") as fh:
                    json.dump(payload, fh)
                os.replace(tmp, path)
            except Exception:
                pass
        return wrapped

    for mod in list(sys.modules.values()):
        cls = getattr(mod, "BatchProgressTracker", None)
        if isinstance(cls, type) and not getattr(cls, "_cellenics_patched", False):
            try:
                cls.complete_task = make_wrapped(cls.complete_task)
                cls._cellenics_patched = True
            except Exception:
                pass

_cellenics_patch_batch_progress()
)"


#' Enable per-cluster CASSIA progress reporting
#'
#' Points CASSIA's batch tracker at `progress_file` (the generic per-task
#' progress file) and installs the runtime hook on the tracker. Best-effort.
#'
#' @param progress_file character path to the shared progress file
#' @export
#'
enable_cassia_batch_progress <- function(progress_file) {
  tryCatch({
    reticulate::py_run_string(sprintf(
      'import os; os.environ["CASSIA_PROGRESS_FILE"] = "%s"',
      progress_file
    ))
    reticulate::py_run_string(CASSIA_BATCH_PROGRESS_PATCH)
  }, error = function(e) invisible(NULL))
}


#' Disable per-cluster CASSIA progress reporting
#'
#' Unsets `CASSIA_PROGRESS_FILE` so the (still-installed) runtime hook becomes a
#' no-op after the run. Best-effort.
#'
#' @export
#'
disable_cassia_batch_progress <- function() {
  tryCatch(
    reticulate::py_run_string(
      'import os; os.environ.pop("CASSIA_PROGRESS_FILE", None)'
    ),
    error = function(e) invisible(NULL)
  )
}


#' Compute per-cluster marker genes for CASSIA
#'
#' Computes one-vs-rest marker genes per cluster via [computeMarkerStats()]
#' (presto::wilcoxauc on a per-cluster downsample, the same engine and
#' downsampling used by [getTopMarkerGenes()]). presto's statistics are then
#' converted to the exact shape and default filtering of a
#' `Seurat::FindAllMarkers(only.pos = TRUE)` call (see
#' [presto_markers_as_find_all_markers()]), which is the format CASSIA expects.
#' This avoids running FindAllMarkers directly, which chokes on the
#' multi-assay / spatial Seurat objects (see git history / DietSeurat notes).
#' When gene annotations are available, gene ids are mapped to symbols so the
#' LLM receives human readable names.
#'
#' @param data Seurat object with the `cluster_col` metadata column
#' @param max_cells_per_cluster integer. Max cells sampled per cluster
#' @param cluster_col string. Metadata column holding the cluster labels to
#'  annotate (the labels become the `cluster` column of the returned markers)
#'
#' @return data.frame with FindAllMarkers columns:
#'  p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
#' @export
#'
get_cassia_markers <- function(
  data,
  max_cells_per_cluster = 1000,
  cluster_col = "seurat_clusters"
) {
  clusters <- data@meta.data[[cluster_col]]
  cluster_levels <- unique(clusters[!is.na(clusters)])

  # one cell set (of cell_ids) per cluster, in cluster_levels order; the group
  # index i returned by computeMarkerStats maps back to cluster_levels[i]
  cell_sets_ids <- lapply(cluster_levels, function(cl) {
    data@meta.data$cells_id[which(clusters == cl)]
  })

  stats <- computeMarkerStats(
    data,
    cell_sets_ids,
    downsample_n = max_cells_per_cluster,
    seed = ULTIMATE_SEED
  )

  markers <- presto_markers_as_find_all_markers(
    stats,
    cluster_labels = cluster_levels,
    n_total_genes = nrow(data[["RNA"]])
  )
  message("CASSIAAnnotate: found ", nrow(markers), " markers across clusters")

  markers <- add_marker_gene_symbols(markers, data)

  return(markers)
}


#' Convert presto::wilcoxauc stats to FindAllMarkers output
#'
#' Reproduces the output of `Seurat::FindAllMarkers(only.pos = TRUE)` with
#' default parameters from presto statistics computed on the same cells.
#'
#' The wilcoxon p-value, `pct.1` (= `pct_in / 100`) and `pct.2`
#' (= `pct_out / 100`) come straight from presto (FindAllMarkers uses presto
#' internally, so these match). `avg_log2FC` is NOT taken from presto's `logFC`
#' (a natural-log difference of means-of-logs); it is recomputed the Seurat way
#' via [compute_seurat_log2fc()] (log2 of the ratio of mean expm1 expression,
#' pseudocount 1). The default FindAllMarkers filters are then applied:
#' `max(pct.1, pct.2) >= min_pct`, `abs(avg_log2FC) >= logfc_threshold`,
#' `p_val < return_thresh`, and (for only.pos) `avg_log2FC > 0`. `p_val_adj` is
#' Bonferroni across all genes in the object.
#'
#' @param stats list from [computeMarkerStats()] (markers, mat, groups)
#' @param cluster_labels character vector mapping group index -> cluster name
#' @param n_total_genes int total genes in the object (Bonferroni denominator)
#' @param logfc_threshold,min_pct,return_thresh FindAllMarkers defaults
#'
#' @return data.frame with columns
#'  p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
#' @export
#'
presto_markers_as_find_all_markers <- function(
  stats,
  cluster_labels,
  n_total_genes,
  logfc_threshold = 0.1,
  min_pct = 0.01,
  return_thresh = 0.01
) {
  all_markers <- stats$markers

  # Seurat-style log2 fold change on the same downsampled cells
  log2fc <- compute_seurat_log2fc(stats$mat, stats$groups)
  all_markers$avg_log2FC <- log2fc[
    cbind(all_markers$feature, as.character(all_markers$group))
  ]

  markers <- all_markers |>
    dplyr::transmute(
      p_val = pval,
      avg_log2FC = avg_log2FC,
      pct.1 = pct_in / 100,
      pct.2 = pct_out / 100,
      p_val_adj = pmin(pval * n_total_genes, 1),
      cluster = cluster_labels[group],
      gene = feature
    ) |>
    # default FindAllMarkers(only.pos = TRUE) filtering
    dplyr::filter(
      pmax(pct.1, pct.2) >= min_pct,
      abs(avg_log2FC) >= logfc_threshold,
      p_val < return_thresh,
      avg_log2FC > 0
    ) |>
    dplyr::arrange(cluster, p_val, dplyr::desc(avg_log2FC))

  rownames(markers) <- NULL
  return(markers)
}


#' Seurat-style log2 fold change (one-vs-rest)
#'
#' Reproduces `Seurat::FoldChange()` for log-normalized data (its
#' `log1pdata.mean.fxn`): for each group,
#' `log2((sum(expm1(x_in)) + pc) / n_in) - log2((sum(expm1(x_out)) + pc) / n_out)`.
#' Note the pseudocount `pc` (default 1) is added to the summed expression and
#' then divided by the cell count, i.e. it is `mean + pc/n`, NOT `mean + pc`.
#'
#' @param mat dgCMatrix, genes x cells (log-normalized `data` layer)
#' @param groups vector of group labels, one per column of `mat`
#' @param pseudocount numeric pseudocount added to the summed expression
#' @param base numeric log base (2 for avg_log2FC)
#'
#' @return numeric matrix genes x groups of log2 fold changes; rownames are
#'  gene ids, colnames are the (as.character) group labels
#' @export
#'
compute_seurat_log2fc <- function(mat, groups, pseudocount = 1, base = 2) {
  # expm1 on nonzero entries only (expm1(0) == 0 keeps the matrix sparse)
  expm1_mat <- mat
  expm1_mat@x <- expm1(expm1_mat@x)

  total_sum <- Matrix::rowSums(expm1_mat)
  n_total <- ncol(expm1_mat)
  group_levels <- sort(unique(groups))

  fc <- vapply(group_levels, function(g) {
    in_cols <- which(groups == g)
    sum_in <- Matrix::rowSums(expm1_mat[, in_cols, drop = FALSE])
    n_in <- length(in_cols)
    sum_out <- total_sum - sum_in
    n_out <- n_total - n_in
    log((sum_in + pseudocount) / n_in, base) -
      log((sum_out + pseudocount) / n_out, base)
  }, numeric(nrow(expm1_mat)))

  rownames(fc) <- rownames(mat)
  colnames(fc) <- as.character(group_levels)
  return(fc)
}


#' Map marker gene ids to gene symbols
#'
#' CASSIA needs gene symbols to reason about cell types. If the Seurat object
#' carries `gene_annotations` (created in the pipeline), map the marker `gene`
#' column from feature id to `original_name`. If no annotations are present 
#' (or a feature has no symbol) the original value is kept unchanged.
#'
#' @param markers data.frame of markers with a `gene` column (FindAllMarkers
#'  shaped; see [get_cassia_markers()])
#' @param data Seurat object
#'
#' @return markers data.frame with the `gene` column mapped to symbols
#' @export
#'
add_marker_gene_symbols <- function(markers, data) {
  annot <- data@misc$gene_annotations
  if (is.null(annot) || !all(c("input", "original_name") %in% colnames(annot))) {
    return(markers)
  }

  symbol_by_id <- stats::setNames(annot$original_name, annot$input)
  mapped <- symbol_by_id[markers$gene]
  markers$gene <- ifelse(is.na(mapped), markers$gene, unname(mapped))

  return(markers)
}


#' Run the CASSIA annotation pipeline
#'
#' Runs [CASSIA::runCASSIA_pipeline()] on the supplied markers using a single
#' OpenRouter model for every agent, then merges the cluster-level predictions
#' back onto the Seurat metadata with [CASSIA::add_cassia_to_seurat()]. Each of
#' the (up to three) progressively coarser merged groupings is stored in a
#' `CASSIA_merged_grouping_1` / `_2` / `_3` column.
#' Modelled after `data-raw/testing_cassia.R`.
#'
#' @param data Seurat object with the `cluster_col` metadata column
#' @param markers data.frame of markers from [get_cassia_markers()]
#' @param tissue string. Free text tissue type
#' @param species string. Free text species
#' @param model string. Bedrock model id used for every CASSIA agent
#' @param provider string. CASSIA provider — the OpenAI-compatible Bedrock
#'  endpoint URL (see [set_cassia_bedrock_auth()])
#' @param additional_info string or NULL. Optional free-text study context passed
#'  to the CASSIA agents; blank/whitespace is treated as no context
#' @param cluster_col string. Metadata column whose labels match the marker
#'  `cluster` column, used to map CASSIA predictions back onto cells
#' @param progress_file character. When set, CASSIA's per-cluster batch progress
#'  is mirrored to this shared progress file so the UI can show fine-grained
#'  progress as clusters complete
#'
#' @return Seurat object with CASSIA predictions in the metadata slot
#' @export
#'
run_cassia <- function(
  data, markers, tissue, species, model, provider,
  additional_info = NULL,
  cluster_col = "seurat_clusters",
  progress_file = NULL
) {
  # Treat blank/whitespace context as "none" so CASSIA doesn't inject an empty
  # additional-info line into the prompts.
  if (!is.null(additional_info) && !nzchar(trimws(additional_info))) {
    additional_info <- NULL
  }
  output_dir <- tempfile("cassia_results_")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Surface CASSIA's per-cluster batch progress to the UI as each cluster
  # finishes. CASSIA exposes no progress callback and writes its results only at
  # the end, so we hook its in-memory progress tracker at runtime (see
  # enable_cassia_batch_progress). Best-effort: if the hook can't attach we
  # simply fall back to the coarser phase-level progress.
  if (!is.null(progress_file)) {
    enable_cassia_batch_progress(progress_file)
    on.exit(disable_cassia_batch_progress(), add = TRUE)
  }

  # One worker per cluster so all clusters annotate concurrently, capped at 50
  # to stay within Bedrock's per-model throughput limits on large datasets.
  n_clusters <- length(unique(markers$cluster))
  max_workers <- min(n_clusters, 50)

  tryCatch(
    CASSIA::runCASSIA_pipeline(
      output_file_name = "results",
      tissue = tissue,
      species = species,
      marker = markers,
      additional_info = additional_info,
      max_workers = max_workers,
      overall_provider = provider,
      annotation_model = model,
      annotationboost_model = model,
      score_model = model,
      merge_model = model,
      output_dir = output_dir,
      # The pre-flight check tests a hard-coded non-Bedrock validation model
      # that this endpoint does not serve; auth is already established above.
      validate_api_keys_before_start = FALSE
    ),
    error = function(e) {
      stop(generateErrorMessage(error_codes$CASSIA_ERROR, e$message))
    }
  )

  # the pipeline writes the final table to
  # <output_dir>/CASSIA_Pipeline_*/03_csv_files/*_FINAL_RESULTS.csv
  results_csv <- list.files(
    output_dir,
    pattern = "_FINAL_RESULTS\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(results_csv) == 0) {
    stop(generateErrorMessage(
      error_codes$CASSIA_ERROR,
      "CASSIA did not produce a results file."
    ))
  }

  data <- CASSIA::add_cassia_to_seurat(
    seurat_obj = data,
    cassia_results_path = results_csv[1],
    cluster_col = cluster_col,
    columns_to_include = 1
  )

  return(data)
}


#' Format a CASSIA cell class for patching through the API
#'
#' Builds a new cell class from one CASSIA merged grouping stored in the
#' `CASSIA_merged_grouping_<grouping>` metadata column, one child cell set per
#' predicted cell type. CASSIA produces up to three progressively coarser
#' groupings; the class key is suffixed with the grouping number (e.g.
#' `CASSIA-<tissue>-<species>-2`). Mirrors [format_sctype_cell_sets()].
#'
#' @param data Seurat object with CASSIA annotations in
#'  `CASSIA_merged_grouping_<grouping>`
#' @param species string. Free text species
#' @param tissue string. Free text tissue type
#' @param grouping integer. Which merged grouping to format (1, 2 or 3)
#'
#' @return list, a cellsets cell class ready to patch through the API
#' @export
#'
format_cassia_cell_sets <- function(data, species, tissue, grouping) {
  annotation_column <- paste0("CASSIA_merged_grouping_", grouping)
  cell_class_key <- paste0("CASSIA-", grouping)

  cell_class <-
    list(
      key = uuid::UUIDgenerate(use.time = TRUE),
      name = cell_class_key,
      rootNode = TRUE,
      type = "cellSets",
      children = list()
    )

  annotations <- data@meta.data[[annotation_column]]
  unique_annotations <- unique(annotations)
  unique_annotations <- unique_annotations[!is.na(unique_annotations)]

  for (annotation in unique_annotations) {
    new_cell_set <- list(
      key = uuid::UUIDgenerate(use.time = TRUE),
      name = annotation,
      rootNode = FALSE,
      type = "cellSets",
      color = sample(data@misc$color_pool, 1),
      cellIds = data@meta.data[which(annotations == annotation), "cells_id"]
    )
    cell_class$children <- append(cell_class$children, list(new_cell_set))
  }

  return(cell_class)
}
