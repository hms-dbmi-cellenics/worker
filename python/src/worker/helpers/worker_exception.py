class WorkerException(Exception):
    def __init__(self, error_code, user_message):
        self.error_code = error_code
        self.user_message = user_message
