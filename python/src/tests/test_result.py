import json

import pytest
from worker.result import Result


class TestResult:
    def test_throws_typeerror_on_empty_result(self):
        with pytest.raises(TypeError):
            Result()

    def test_can_create_result_with_no_type_or_encoding(self):
        Result("")

    def test_result_works_with_default_type_and_encoding(self):
        result = Result("")

        assert result.content_encoding == "utf-8"
        assert result.content_type == "application/json"

    def test_result_works_with_custom_type_and_encoding(self):
        result = Result(
            "", content_type="image/svg+xml", content_encoding="base64"
        )

        assert result.content_encoding == "base64"
        assert result.content_type == "image/svg+xml"

    def test_get_result_object_returns_a_dict(self):
        result = Result("")
        o = result.get_result_object()

        assert type(o) == dict

    def test_get_result_object_returns_custom_type_and_encoding(self):
        result = Result(
            "", content_type="image/svg+xml", content_encoding="base64"
        )

        o = result.get_result_object()

        assert o["content-type"] == "image/svg+xml"
        assert o["content-encoding"] == "base64"

    def test_get_result_object_returns_type_when_resp_format_is_requested(
        self,
    ):
        result = Result("")
        o = result.get_result_object(resp_format=True)

        assert o["type"]

    def test_get_result_object_returns_appropriate_path_when_s3_path_is_patched(
        self,
    ):
        result = Result("")
        o = result.get_result_object(resp_format=True, s3_path="my/s3/path")

        assert o["type"] == "s3-path"
        assert o["body"] == "my/s3/path"

    def test_get_result_length_is_no_smaller_than_length_of_json_text(self):
        result = Result("")
        length = result.get_result_length()

        assert length >= len(json.dumps(result.get_result_object()))
