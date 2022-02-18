class RWorkerException(Exception):
    def __init__(self, user_message, code):
        self.user_message = user_message
        self.code = code
