import json


class Result:
    def __init__(
        self,
        result,
        content_type="application/json",
        content_encoding="utf-8",
        error=False,
        cacheable=True,
    ):
        self.result = result
        self.content_type = content_type
        self.content_encoding = content_encoding
        self.error = error
        self.cacheable = cacheable

    def get_result_object(self, resp_format=False, s3_path=None):
        obj = {
            "content-type": self.content_type,
            "content-encoding": self.content_encoding,
            "body": self.result,
        }

        if resp_format and s3_path:
            obj["type"] = "s3-path"
            obj["body"] = s3_path

        if resp_format and not s3_path:
            obj["type"] = "inline"

        return obj

    def get_result_length(self):
        return len(json.dumps(self.get_result_object()))
