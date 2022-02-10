import json

import botocore.session
import mock
import pytest
from botocore.stub import Stubber

from worker.config import config
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
        stubbed_client = botocore.session.get_session().create_client(
            "s3", **config.BOTO_RESOURCE_KWARGS
        )
        stubber = Stubber(stubbed_client)
        stubber.activate()

        r = Result({})
        resp = Response(self.request, r)
        key = resp._upload(r)

        assert key == self.request["ETag"]

    def test_construct_response_msg_works(self):
        resp = Response(self.request, Result({"result1key": "result1val"}))
        response_msg = resp._construct_response_msg()

        assert response_msg["request"] == self.request
        assert response_msg["response"]["cacheable"] == True
        assert response_msg["response"]["error"] == False

    # These tests only work locally with inframock running. Keeping in case
    # we want to mock to be able to run these tests.
    @mock.patch("boto3.client")
    def test_publishing_long_responses_get_pushed_to_s3(self, mocked_client, mocker):
        stubbed_client = botocore.session.get_session().create_client(
            "s3", **config.BOTO_RESOURCE_KWARGS
        )
        stubber = Stubber(stubbed_client)
        stubber.activate()

        result = Result(
            "a" * 512 * 1024,
            content_encoding="base64",
            content_type="application/octet-stream",
        )

        resp = Response(self.request, result)
        spy = mocker.spy(resp, "_upload")

        with mock.patch("worker.response.Emitter") as redis_emitter:
            resp.publish()
            assert redis_emitter.call_count == 1
        assert spy.call_count == 1

    @mock.patch("boto3.client")
    def test_publishing_one_long_response_results_in_both_being_pushed_to_s3(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client(
            "s3", **config.BOTO_RESOURCE_KWARGS
        )
        stubber = Stubber(stubbed_client)
        stubber.activate()

        result = Result(
            "a" * 512 * 1024,
            content_encoding="base64",
            content_type="application/octet-stream",
        )

        resp = Response(self.request, result)
        spy = mocker.spy(resp, "_upload")

        with mock.patch("worker.response.Emitter") as redis_emitter:
            resp.publish()
            assert redis_emitter.call_count == 1
        assert spy.call_count == 1

    @mock.patch("boto3.client")
    def test_old_requests_do_get_sent(self, mocked_client, mocker):
        stubbed_client = botocore.session.get_session().create_client(
            "s3", **config.BOTO_RESOURCE_KWARGS
        )
        stubber = Stubber(stubbed_client)
        stubber.activate()

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
            assert redis_emitter.call_count == 1
        assert spy.call_count == 1
