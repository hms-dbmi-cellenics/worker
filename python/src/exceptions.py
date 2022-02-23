from enum import Enum


class WorkerException(Exception):
    def __init__(self, error_code, user_message):
        self.error_code = error_code
        self.user_message = user_message


class RWorkerException(WorkerException):
    def __init__(self, error_code, user_message):
        super().__init__(error_code, user_message)


class PythonWorkerException(WorkerException):
    def __init__(self, error_code, user_message):
        super().__init__(error_code, user_message)


class ErrorCodes(Enum):
    R_WORKER_ERROR = "R_WORKER_ERROR"
    PYTHON_WORKER_ERROR = "PYTHON_WORKER_ERROR"
    INVALID_INPUT = "INVALID_INPUT"


def raise_if_error(result):
    error = result.get("error")
    if error:
        user_message = error.get("user_message", "Unexpected error occured")
        err_code = error.get("error_code", ErrorCodes.R_WORKER_ERROR)
        raise RWorkerException(err_code, user_message)
