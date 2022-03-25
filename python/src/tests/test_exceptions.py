import pytest
from exceptions import RWorkerException, raise_if_error


class TestRaiseIfError:
    def test_should_not_raise_error_if_there_is_no_error(self):
        result = {"data": "no error"}
        raise_if_error(result)

    def test_should_raise_error_if_there_is_an_error(self):
        result = {
            "error": {"user_message": "Some random error", "error_code": "R_MOCK_ERROR"}
        }

        with pytest.raises(RWorkerException):
            raise_if_error(result)
