import gzip
import io
import json
from unittest import TestCase

import boto3
import mock
from botocore.stub import Stubber
from tests.data.embedding import mock_embedding
from worker.config import config
from worker.helpers.s3 import get_embedding

mock_embedding_etag = "mockEmbeddingETag"

class TestS3:
    def get_s3_stub(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        response = {
            "ContentLength": 10,
            "ContentType": "utf-8",
            "ResponseMetadata": {
                "Bucket": config.RESULTS_BUCKET,
            },
        }

        expected_params = {
            "Bucket": config.RESULTS_BUCKET,
            "Key": mock_embedding_etag,
        }
        stubber = Stubber(s3)
        stubber.add_response("head_object", response, expected_params)

        # Get object
        content_string = json.dumps(mock_embedding).encode("utf-8")
        content_bytes = gzip.compress(content_string)
        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        response = {
            "ContentLength": len(content_bytes),
            "ContentType": "application/gzip",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.RESULTS_BUCKET,
            },
        }
        stubber.add_response("get_object", response, expected_params)
        return (stubber, s3)

    def test_get_embedding_should_not_replace_nulls_if_not_formatted_for_r(self):
      stubber, s3 = self.get_s3_stub()

      na_positions = []
      for idx, val in enumerate(mock_embedding):
        if val is None:
          na_positions.append(idx)

      with mock.patch("boto3.client") as n, stubber:
          n.return_value = s3

          request = get_embedding(mock_embedding_etag, format_for_r=False)

          for idx, val in enumerate(request):
            if idx in na_positions:
              assert val is None
            else:
              assert val is not None

    def test_get_embedding_should_replace_nulls_if_formatted_for_r(self):
      stubber, s3 = self.get_s3_stub()

      na_positions = []
      for idx, val in enumerate(mock_embedding):
        if val is None:
          na_positions.append(idx)

      with mock.patch("boto3.client") as n, stubber:
          n.return_value = s3

          request = get_embedding(mock_embedding_etag, format_for_r=True)

          for idx, val in enumerate(request):
            if idx in na_positions:
              assert val == ['NA', 'NA']
            else:
              assert val is not None