import base64
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
            "signedUrl": "mockSignedUrl",
            "requestProps": {}
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

    def test_construct_response_msg_works_with_signed_url(self):
        resp = Response(self.request, Result({"result1key": "result1val"}))
        response_msg = resp._construct_response_msg()

        assert response_msg["request"] == self.request
        assert response_msg["response"]["cacheable"] is True
        assert response_msg["response"]["error"] is False
        assert response_msg["response"]["signedUrl"] is "mockSignedUrl"

    def test_construct_response_msg_works_with_data(self):
        resp = Response(self.request, Result({"result1key": "result1val"}))

        data = bytes([1,2,3,4,5])

        response_msg = resp._construct_response_msg(data)

        assert response_msg == base64.b64encode(data)

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

    @mock.patch("boto3.client")
    def test_construct_data_for_upload_handles_bytes_data(self, mocked_client):
        """Test that bytes data from R response is handled correctly."""
        # create a result with bytes data (as would come from R)
        bytes_data = b"binary data from R response"
        result = Result(bytes_data)
        result.data = bytes_data

        resp = Response(self.request, result)

        # mock the Emitter to prevent Redis connection
        with mock.patch("worker.response.Emitter"):
            # should not raise an error
            compressed_body, compressed_bytes = resp._construct_data_for_upload()

            # verify compression was performed
            assert compressed_body is not None
            assert compressed_bytes is not None
            assert len(compressed_bytes) > 0

    @mock.patch("boto3.client")
    def test_construct_data_for_upload_prefers_bytes_over_string(
        self, mocked_client
    ):
        """Test that bytes data takes precedence over string data."""
        # create result with bytes
        bytes_data = b"bytes data"
        result = Result("")
        result.data = bytes_data

        resp = Response(self.request, result)

        # mock the compression to track what data is compressed
        compressed_data = None

        def mock_compress(data):
            nonlocal compressed_data
            compressed_data = data
            import io
            import gzip
            buf = io.BytesIO()
            with gzip.GzipFile(fileobj=buf, mode="wb") as f:
                f.write(data if isinstance(data, bytes) else data.encode())
            return buf

        # mock the Emitter to prevent Redis connection
        with mock.patch("worker.response.Emitter"):
            with mock.patch("gzip.compress", side_effect=mock_compress):
                compressed_body, compressed_bytes = (
                    resp._construct_data_for_upload()
                )

            # verify bytes data was used (not string)
            assert compressed_data == bytes_data

    @mock.patch("boto3.client")
    def test_construct_data_for_upload_handles_string_data(self, mocked_client):
        """Test that string data is correctly encoded to bytes."""
        string_data = "text data from result"
        result = Result(string_data)

        resp = Response(self.request, result)

        # mock the Emitter to prevent Redis connection
        with mock.patch("worker.response.Emitter"):
            compressed_body, compressed_bytes = resp._construct_data_for_upload()

            # verify compression was performed
            assert compressed_body is not None
            assert compressed_bytes is not None
            assert len(compressed_bytes) > 0

    @mock.patch("boto3.client")
    def test_construct_data_for_upload_logs_compression_stats(
        self, mocked_client, mocker
    ):
        """Test that compression stats are logged correctly."""
        data = "x" * 10000  # 10KB of data
        result = Result(data)

        resp = Response(self.request, result)

        # mock the logger to verify it was called
        mock_logger = mocker.patch("worker.response.info")

        # mock the Emitter to prevent Redis connection
        with mock.patch("worker.response.Emitter"):
            compressed_body, compressed_bytes = resp._construct_data_for_upload()

            # verify that compression statistics were logged
            # should have at least one call about compression
            assert mock_logger.call_count >= 2
            # find the call with compression size info
            log_calls = [str(call) for call in mock_logger.call_args_list]
            assert any("Compressed from" in str(call) for call in log_calls)
