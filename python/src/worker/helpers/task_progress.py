import json
import os
from concurrent.futures import (
    ThreadPoolExecutor,
    TimeoutError as FuturesTimeoutError,
)

from socket_io_emitter import Emitter
from worker_status_codes import STARTED_TASK

from ..config import config

# how often (seconds) to poll the progress file and emit a heartbeat
PROGRESS_POLL_SECONDS = 5


def progress_file_path(experiment_id, task_name):
    """Conventional path for a task's progress file on the shared /data volume.

    Keyed by task name (the work request name, which equals the R handler /
    python task class name), so this is not specific to any single task.
    """
    return f"/data/{experiment_id}/{task_name}.progress.json"


def _read_progress(path, default):
    # Best-effort: read {step, total, message} written by the R worker. Returns
    # (message, percent) where percent is an int 0-100 when step/total are
    # present (so the UI can show a progress bar), otherwise None. Falls back to
    # `default` with no percent.
    try:
        with open(path) as f:
            progress = json.load(f)
    except (OSError, ValueError):
        return default, None

    step, total = progress.get("step"), progress.get("total")
    text = progress.get("message") or default
    if step is not None and total:
        pct = round(100 * step / total)
        return f"{text} ({step}/{total})", pct
    return text, None


def run_with_progress(
    experiment_id,
    task_name,
    target,
    default_message="Working...",
    poll_seconds=PROGRESS_POLL_SECONDS,
):
    """Run `target` (a no-arg callable) while forwarding progress to the UI.

    `target` typically blocks for a long time (e.g. an HTTP call to the R
    worker). It is run in a background thread; meanwhile we poll the task's
    progress file every `poll_seconds` and emit a heartbeat with the current
    status so the UI can show intermediate progress. Returns `target`'s result
    (re-raising any exception it raised).
    """
    io = Emitter({"client": config.REDIS_CLIENT})
    path = progress_file_path(experiment_id, task_name)

    # clear any stale progress from a previous run
    try:
        os.remove(path)
    except OSError:
        pass

    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(target)
            while True:
                try:
                    return future.result(timeout=poll_seconds)
                except FuturesTimeoutError:
                    message, percent = _read_progress(path, default_message)
                    io.Emit(
                        f"Heartbeat-{experiment_id}",
                        {
                            "type": "WorkResponse",
                            "status_code": STARTED_TASK,
                            "user_message": message,
                            # int 0-100 for a progress bar, or null for
                            # message-only phases (resets any prior bar)
                            "progress_percent": percent,
                        },
                    )
    finally:
        try:
            os.remove(path)
        except OSError:
            pass
