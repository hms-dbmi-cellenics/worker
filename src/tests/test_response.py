from response import Response
from result import Result
import boto3
import botocore.session
from botocore.stub import Stubber
import pytest
from config import get_config
import json
import mock

config = get_config()


class TestResponse:
    def test_throws_on_empty_response_init(self):
        with pytest.raises(TypeError):
            Response()

    def test_throws_on_missing_results(self):
        with pytest.raises(TypeError):
            Response({})

    @mock.patch("boto3.client")
    def test_upload_returns_key_with_uuid_as_folder_when_uploading(self, mocked_client):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        r = Result("{ }")
        resp = Response(request, [r])
        key = resp._upload(r)
        key_folder, *rest = key.split("/")

        assert key_folder == request["uuid"]

    @mock.patch("boto3.client")
    def test_send_notification_pushes_notification_to_sns(self, mocked_client):
        response = {
            "MessageId": "83a8d61a-1056-5e9c-972e-8134130e1d1b",
            "ResponseMetadata": {
                "RequestId": "ac541ad0-560f-5cb4-8382-4dfe5557ff33",
                "HTTPStatusCode": 200,
                "HTTPHeaders": {
                    "x-amzn-requestid": "ac541ad0-560f-5cb4-8382-4dfe5557ff33",
                    "content-type": "text/xml",
                    "content-length": "294",
                    "date": "Thu, 07 May 2020 12:37:43 GMT",
                },
                "RetryAttempts": 0,
            },
        }

        result_object = {"obj_key": "obj_value"}

        request = {
            "TargetArn": "arn:aws:sns:{}:{}:{}".format(
                config.AWS_REGION, config.AWS_ACCOUNT_ID, config.get_sns_topic()
            ),
            "Message": json.dumps({"default": json.dumps(result_object)}),
            "MessageStructure": "json",
        }

        stubbed_client = botocore.session.get_session().create_client("sns")
        stubber = Stubber(stubbed_client)
        stubber.add_response("publish", response, request)
        stubber.activate()
        mocked_client.return_value = stubbed_client

        request = {"uuid": "totally-random-uuid"}
        resp = Response(request, [])

        resp._send_notification(result_object)
        stubber.assert_no_pending_responses()

    def test_get_response_msg_returns_original_request_object(self):
        request = {"uuid": "totally-random-uuid"}
        results = [
            Result({"result1key": "result1val"}),
            Result({"result2key": "result2val"}),
        ]
        resp = Response(request, results)

        print(resp._get_response_msg())

        assert resp._get_response_msg()["request"] == request

    def test_get_response_msg_returns_result_object_definitions(self):
        request = {"uuid": "totally-random-uuid"}
        results = [
            Result({"result1key": "result1val"}, content_encoding="base64"),
            Result({"result2key": "result2val"}, content_encoding="base64"),
        ]

        resp = Response(request, results)
        results_msg = resp._get_response_msg()["results"]

        for msg in results_msg:
            assert msg["content-encoding"] == "base64"
            assert msg["type"] == "inline"

    @mock.patch("boto3.client")
    def test_publishing_long_responses_get_pushed_to_s3(self, mocked_client, mocker):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a" * 512 * 1024,
                content_encoding="base64",
                content_type="application/octet-stream",
            ),
            Result(
                "a" * 512 * 1024,
                content_encoding="base64",
                content_type="application/octet-stream",
            ),
        ]

        resp = Response(request, results)
        spy = mocker.spy(resp, "_upload")

        resp.publish()
        assert spy.call_count == 2

    @mock.patch("boto3.client")
    def test_publishing_one_long_response_results_in_both_being_pushed_to_s3(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a" * 512 * 1024,
                content_encoding="base64",
                content_type="application/octet-stream",
            ),
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
        ]

        resp = Response(request, results)
        spy = mocker.spy(resp, "_upload")

        resp.publish()
        assert spy.call_count == 2

    @mock.patch("boto3.client")
    def test_publishing_one_short_file_results_in_no_s3_uploads(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
        ]

        resp = Response(request, results)
        spy = mocker.spy(resp, "_upload")

        resp.publish()
        assert spy.call_count == 0

    @mock.patch("boto3.client")
    def test_publishing_multiple_short_files_results_in_no_s3_uploads(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
        ]

        resp = Response(request, results)
        spy = mocker.spy(resp, "_upload")

        resp.publish()
        assert spy.call_count == 0

    @mock.patch("boto3.client")
    def test_publishing_a_short_file_results_in_inlined_response(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
            Result(
                "a", content_encoding="base64", content_type="application/octet-stream"
            ),
        ]

        resp = Response(request, results)
        result = resp.publish()

        assert "inline" in result
        assert "s3-path" not in result

    @mock.patch("boto3.client")
    def test_publishing_a_long_file_reuslts_in_s3_path_response(
        self, mocked_client, mocker
    ):
        stubbed_client = botocore.session.get_session().create_client("s3")
        stubber = Stubber(stubbed_client)
        stubber.activate()

        request = {"uuid": "totally-random-uuid"}
        results = [
            Result(
                "a" * 512 * 1024,
                content_encoding="base64",
                content_type="application/octet-stream",
            ),
        ]

        resp = Response(request, results)
        result = resp.publish()

        assert "inline" not in result
        assert "s3-path" in result
