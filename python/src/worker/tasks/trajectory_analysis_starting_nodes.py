import json

import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.s3 import get_embedding
from ..result import Result
from . import Task


class GetTrajectoryAnalysisStartingNodes(Task):
    def __init__(self, msg):
        super().__init__(msg)

    def _format_result(self, result):
        return Result(result["data"])

    def _format_request(self):

        embedding_etag = self.task_def["embedding"]["ETag"]
        embedding = get_embedding(embedding_etag, format_for_r=True)

        request = {
            "embedding": embedding,
            "embedding_settings": {
                "method": self.task_def["embedding"]["method"],
                "methodSettings": self.task_def["embedding"]["methodSettings"]
            },
            "clustering_settings": {
                "method": self.task_def["clustering"]["method"],
                "resolution": self.task_def["clustering"]["resolution"],
            },
        }

        return request

    @xray_recorder.capture("GetTrajectoryAnalysisStartingNodes.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/runTrajectoryAnalysisStartingNodesTask",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError occurred. otherwise r.json() will fail
        r.raise_for_status()
        # The index order relies on cells_id in an ascending form. The order is made in the R part.
        result = r.json()
        raise_if_error(result)

        return self._format_result(result)
