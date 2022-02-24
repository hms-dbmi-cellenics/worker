import traceback
from logging import error

from aws_xray_sdk.core import xray_recorder


def xray_log_exception(task, error_object):
    trace = (
        f"Exception for task {task.__class__.__name__}:\n" f"{traceback.format_exc()}"
    )

    error(trace)
    xray_recorder.current_segment().add_exception(error_object, trace)
