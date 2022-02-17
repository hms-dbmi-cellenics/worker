class RWorkerException(Exception):
    def __init__(self, message, code):
        self.message = message
        self.code = code
