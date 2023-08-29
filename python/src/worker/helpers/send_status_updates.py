from worker_status_codes import (
    DOWNLOAD_EXPERIMENT,
    LOAD_EXPERIMENT,
    STARTED_TASK,
    COMPRESSING_TASK_DATA,
    UPLOADING_TASK_DATA,
    FINISHED_TASK,
)


# format_task_name is a helper function that takes a request and returns a formatted task name
# if the request has a name, otherwise it returns an empty string
def format_task_name(request):
    try:
        task_name = request["body"]["name"]
    except TypeError:
        return ""
    # Remove "Get" from the name
    name_without_get = task_name.replace("Get", "")

    # Insert a space before capital letters, then strip leading/trailing spaces
    final_name = "".join(
        [c if not c.isupper() else f" {c}" for c in name_without_get]
    ).strip()

    return final_name


def format_user_message(status_code, request):
    task_name = format_task_name(request)
    if status_code == DOWNLOAD_EXPERIMENT or status_code == LOAD_EXPERIMENT:
        return "Accessing the Seurat object for your analysis"
    if status_code == STARTED_TASK:
        return f"Working on the requested task {task_name}"
    if status_code == COMPRESSING_TASK_DATA:
        return f"Compressing the results for the requested task {task_name}"
    if status_code == UPLOADING_TASK_DATA:
        return f"Finalizing results for the requested task {task_name}"
    if status_code == FINISHED_TASK:
        return f"Displaying results for the requested task {task_name}"
    return status_code


def send_status_update(io, experiment_id, status_code, request=None):
    io.Emit(
        f"Heartbeat-{experiment_id}",
        {
            "type": "WorkResponse",
            "status_code": status_code,
            "user_message": format_user_message(status_code, request),
        },
    )
