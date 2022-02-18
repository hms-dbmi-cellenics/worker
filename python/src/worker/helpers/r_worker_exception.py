class RWorkerException(Exception):
    def __init__(self, user_message, error_code):
        self.user_message = user_message
        self.error_code = error_code
