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