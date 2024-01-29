import io
import json
import os
import gzip

import boto3
import mock
import pandas as pd
import pytest
import responses
from botocore.stub import Stubber
from exceptions import RWorkerException
from tests.data.cell_set_types import cell_set_types
from worker.config import config
from tests.data.embedding import mock_embedding
from worker.tasks.download_annot_seurat_object import DownloadAnnotSeuratObject

import pdb
from worker.helpers.s3 import get_cell_sets, get_embedding
from worker.helpers.cell_sets_dict import get_cell_sets_dict_for_r

mock_embedding_etag = "mockEmbeddingETag"

class TestDownloadAnnotSeuratObject:
    @pytest.fixture(autouse=True)
    def get_request(self):
        self.correct_request = {
            "experimentId": "e52b39624588791a7889e39c617f669e",
            "timeout": "2099-12-31 00:00:00",
            "Authorization" : "mock_authJwt",
            "body": {
                "name": "DownloadAnnotSeuratObject",
                "embeddingETag": mock_embedding_etag,
                "embeddingMethod": "umap",
                "isSeurat": False
            }
        }

    def get_s3_stub(self):
        s3 = boto3.client("s3", **config.BOTO_RESOURCE_KWARGS)
        stubber = Stubber(s3)

        # Stubbing responses for cell sets head object
        cell_sets_head_object = {
            "params": {
              "Bucket": config.CELL_SETS_BUCKET,
              "Key": config.EXPERIMENT_ID,
            },
            "response": {
              "ContentLength": 10,
              "ContentType": "utf-8",
              "ResponseMetadata": {
                  "Bucket": config.CELL_SETS_BUCKET,
              },
            }
        }

        stubber.add_response("head_object", cell_sets_head_object["response"], cell_sets_head_object["params"])

        # Stubbing responses for cell sets get object
        data = io.BytesIO()
        with open(os.path.join("tests/data", "MockCellSet.json"), "rb") as f:
          content_bytes = f.read()

        data.write(content_bytes)
        data.seek(0)

        cell_sets_get_object = {
            "params": {
              "Bucket": config.CELL_SETS_BUCKET,
              "Key": config.EXPERIMENT_ID,
            },
            "response": {
              "ContentLength": len(content_bytes),
              "ContentType": "utf-8",
              "Body": data,
              "ResponseMetadata": {
                  "Bucket": config.CELL_SETS_BUCKET,
              },
            }
        }

        stubber.add_response("get_object", cell_sets_get_object["response"], cell_sets_get_object["params"])

        embedding_head_object = {
          "params": {
            "Bucket": config.RESULTS_BUCKET,
            "Key": mock_embedding_etag,
          },
          "response": {
              "ContentLength": 10,
              "ContentType": "utf-8",
              "ResponseMetadata": {
                  "Bucket": config.RESULTS_BUCKET,
              },
          }
        }

        # Stubbing response for embedding head object
        stubber.add_response("head_object", embedding_head_object["response"], embedding_head_object["params"])

        # Stubbing response for embedding get object
        content_string = json.dumps(mock_embedding).encode("utf-8")
        content_bytes = gzip.compress(content_string)
        data = io.BytesIO()
        data.write(content_bytes)
        data.seek(0)

        embedding_get_object = {
          "params": {
            "Bucket": config.RESULTS_BUCKET,
            "Key": mock_embedding_etag,
          },
          "response": {
            "ContentLength": len(content_bytes),
            "ContentType": "application/gzip",
            "Body": data,
            "ResponseMetadata": {
                "Bucket": config.RESULTS_BUCKET,
            },
          }
        }

        stubber.add_response("get_object", embedding_get_object["response"], embedding_get_object["params"])

        return (stubber, s3)

    def test_works_with_request(self):
        DownloadAnnotSeuratObject(self.correct_request)

    def test_throws_on_missing_parameters(self):
        with pytest.raises(TypeError):
            DownloadAnnotSeuratObject()

    def test_generates_correct_request_keys(self):
        stubber, s3 = self.get_s3_stub()

        with mock.patch("boto3.client") as n, stubber:
            n.return_value = s3
            request = DownloadAnnotSeuratObject(self.correct_request)._format_request()
            assert isinstance(request, dict)

            # all expected keys are in the request
            expected_keys = [
                "embedding",
                "cellSets",
                "embeddingMethod",
                ]
            assert all(key in request for key in expected_keys)
            stubber.assert_no_pending_responses()
 
    @responses.activate
    def test_should_throw_exception_on_r_worker_error(self):

         error_code = "MOCK_R_WORKER_ERROR"
         user_message = "Some worker error"

         stubber, s3 = self.get_s3_stub()

         payload = {"error": {"error_code": error_code, "user_message": user_message}}

         responses.add(
             responses.POST,
             f"{config.R_WORKER_URL}/v0/DownloadAnnotSeuratObject",
             json=payload,
             status=200,
         )

         with mock.patch("boto3.client") as n, stubber:
             n.return_value = s3
             with pytest.raises(RWorkerException) as exception_info:
                 DownloadAnnotSeuratObject(self.correct_request).compute()

             assert exception_info.value.args[0] == error_code
             assert exception_info.value.args[1] == user_message