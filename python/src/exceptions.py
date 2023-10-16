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


def raise_if_error(result):
    """Raise exception if result is an error.

    Args:
        result (Dict): JSON result returned from the R Worker

    Raises:
        RWorkerException: Raise R worker exception containing the error code
            and user message of the error
    """
    error = result.get("error")
    if error:
        raise RWorkerException(error["error_code"], error["user_message"])
