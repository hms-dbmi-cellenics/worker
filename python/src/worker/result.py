class Result:
    def __init__(
        self,
        data,
        content_type="application/json",
        content_encoding="utf-8",
        error=False,
        cacheable=True,
        upload=True,
    ):
        self.data = data
        self.content_type = content_type
        self.content_encoding = content_encoding
        self.error = error
        self.cacheable = cacheable
        self.upload = upload
