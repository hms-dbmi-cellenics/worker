error_codes <- list(
  GENE_NOT_FOUND = "R_WORKER_GENE_NOT_FOUND",
  EMPTY_CELL_SET = "R_WORKER_EMPTY_CELL_SET",
  COLUMN_NOT_FOUND = "R_WORKER_COLUMN_NOT_FOUND",
  NO_MARKER_GENES = "R_WORKER_NO_MARKER_GENES",
  EMPTY_ROOT_NODES = "R_WORKER_TRAJECTORY_ANALYSIS_EMPTY_ROOT_NODES",
  NO_GENE_SYMBOLS = "R_WORKER_NO_GENE_SYMBOLS",
  UNHANDLED_ERROR = "R_WORKER_UNHANDLED_ERROR"
)

ULTIMATE_SEED <- 42

QUANTILE_THRESHOLD <- 0.95

# annotation type constants
SYM_IDS <- "sym_ids"
SYM_SYM <- "sym_sym"
IDS_SYM <- "ids_sym"
IDS_IDS <- "ids_ids"

RDS_PATH <- "/data/processed.rds"

INTERNAL_RESULTS_PATH <- "/data/rResult"

usethis::use_data(error_codes,
                  ULTIMATE_SEED,
                  QUANTILE_THRESHOLD,
                  SYM_IDS,
                  SYM_SYM,
                  IDS_SYM,
                  IDS_IDS,
                  RDS_PATH,
                  INTERNAL_RESULTS_PATH,
                  internal = TRUE,
                  overwrite = TRUE)
