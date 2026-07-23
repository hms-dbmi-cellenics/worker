#' Task progress reporting
#'
#' Generic (task-agnostic) helpers for surfacing intermediate progress from a
#' long-running worker task to the UI. A task writes `{step, total, message}` to
#' a conventionally named file on the shared `/data` volume; the python worker
#' polls that file while the R task runs and forwards the progress to the UI as
#' heartbeats (see worker/python/src/worker/helpers/task_progress.py).
#'
#' The file is keyed by task name so it is not specific to any single task. The
#' task name matches the work request `name` (`req$body$name`), which the python
#' side derives from the task class name - these are always equal.


#' Path to a task's progress file
#'
#' @param experiment_id character
#' @param task_name character, the work request name (e.g. "CASSIAAnnotate")
#'
#' @return character path on the shared /data volume
#' @export
#'
worker_progress_file <- function(experiment_id, task_name) {
  file.path("/data", experiment_id, paste0(task_name, ".progress.json"))
}


#' Write task progress
#'
#' Best-effort: any failure is swallowed so progress reporting never breaks the
#' task itself. When `step`/`total` are omitted, only the message is written and
#' the UI shows it without a percentage (use for fast prep phases so it doesn't
#' look e.g. 75% done before the long phase starts). Provide `step`/`total` only
#' for phases where a percentage is meaningful.
#'
#' @param experiment_id character
#' @param task_name character, the work request name
#' @param message character human readable phase description
#' @param step,total integer current step and total; omit for message-only
#' @export
#'
write_worker_progress <- function(
  experiment_id, task_name, message, step = NULL, total = NULL
) {
  payload <- list(message = message)
  if (!is.null(step) && !is.null(total)) {
    payload$step <- step
    payload$total <- total
  }
  tryCatch(
    jsonlite::write_json(
      payload,
      worker_progress_file(experiment_id, task_name),
      auto_unbox = TRUE
    ),
    error = function(e) invisible(NULL)
  )
}


#' Remove a task's progress file (best-effort)
#'
#' @param experiment_id character
#' @param task_name character, the work request name
#' @export
#'
clear_worker_progress <- function(experiment_id, task_name) {
  tryCatch(
    unlink(worker_progress_file(experiment_id, task_name)),
    error = function(e) invisible(NULL)
  )
}
