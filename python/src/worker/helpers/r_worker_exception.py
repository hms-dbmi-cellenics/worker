class RWorkerException(Exception):
    def __init__(self, user_msg, error):
        self.user_msg = user_msg
        self.error = error
