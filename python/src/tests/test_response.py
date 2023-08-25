import mock
import pytest

from worker.response import Response
from worker.result import Result


class TestResponse:
    @pytest.fixture(autouse=True)
    def load_correct_definition(self):
        self.request = {
            "body": {
                "name": "DifferentialExpression",
            },
            "experimentId": "random-experiment-id",
            "timeout": "2099-12-31 00:00:00",
            "uuid": "random-uuid",
            "ETag": "random-etag",
            "socketId": "random-socketId",
        }

    def test_throws_on_empty_response_init(self):
        with pytest.raises(TypeError):
            Response()

    def test_throws_on_missing_results(self):
        with pytest.raises(TypeError):
            Response({})

    @mock.patch("boto3.client")
    def test_upload_returns_etag_as_key_when_uploading(self, mocked_client):
        r = Result({})
        resp = Response(self.request, r)
        type = "obj"
        with mock.patch("worker.response.Emitter") as redis_emitter:
            key = resp._upload(r, type)
            assert redis_emitter.call_count >= 1
            assert key == self.request["ETag"]

    def test_construct_response_msg_works(self):
        resp = Response(self.request, Result({"result1key": "result1val"}))
        response_msg = resp._construct_response_msg()

        assert response_msg["request"] == self.request
        assert response_msg["response"]["cacheable"] is True
        assert response_msg["response"]["error"] is False

    @mock.patch("boto3.client")
    def test_publishing_long_responses_get_pushed_to_s3(self, mocked_client, mocker):
        result = Result(
            "a" * 512 * 1024,
            content_encoding="base64",
            content_type="application/octet-stream",
        )

        resp = Response(self.request, result)
        spy = mocker.spy(resp, "_upload")

        with mock.patch("worker.response.Emitter") as redis_emitter:
            resp.publish()
            assert redis_emitter.call_count >= 1
        assert spy.call_count >= 1

    @mock.patch("boto3.client")
    def test_publishing_one_long_response_results_in_both_being_pushed_to_s3(
        self, mocked_client, mocker
    ):
        result = Result(
            "a" * 512 * 1024,
            content_encoding="base64",
            content_type="application/octet-stream",
        )

        resp = Response(self.request, result)
        spy = mocker.spy(resp, "_upload")

        with mock.patch("worker.response.Emitter") as redis_emitter:
            resp.publish()
            assert redis_emitter.call_count >= 1
        assert spy.call_count >= 1

    @mock.patch("boto3.client")
    def test_old_requests_do_get_sent(self, mocked_client, mocker):
        request = self.request
        request["timeout"] = "2000-01-01 00:00:00"

        result = Result(
            "a",
            content_encoding="base64",
            content_type="application/octet-stream",
        )

        resp = Response(self.request, result)
        spy = mocker.spy(resp, "_upload")

        with mock.patch("worker.response.Emitter") as redis_emitter:
            resp.publish()
            assert redis_emitter.call_count >= 1
        assert spy.call_count >= 1
