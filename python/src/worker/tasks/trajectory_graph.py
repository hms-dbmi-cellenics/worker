import gzip
import json

import backoff
import boto3
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..result import Result
from ..tasks import Task


class GetTrajectoryGraph(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        return Result(result["data"])

    def _format_request(self):

        request = {
            "embedding_settings": (
                self.task_def["configureEmbedding"]["embeddingSettings"]
            ),
            "clustering_settings": (
                self.task_def["configureEmbedding"]["clusteringSettings"]
            ),
        }

        return request

    def _recreate_seurat_clusters(self, clusters):

        cluster_array = []
        for cluster in clusters:

            # We get the cluster number by splitting cluster key (e.g. "louvain-1")
            # and getting the last part of the string
            cluster_id = int(cluster["key"].split("-")[-1])

            for cell_id in cluster["cellIds"]:

                cell_id = int(cell_id)

                # cell_id starts from 0, so we have to handle the possibility
                # that the first cell id is a 0
                if len(cluster_array) < cell_id + 1:
                    extend_by = cell_id - len(cluster_array) + 1
                    cluster_array.extend([None] * extend_by)

                cluster_array[cell_id] = cluster_id

        cluster = [cell for cell in cluster_array if cell is not None]
        return cluster

    @xray_recorder.capture("GetTrajectoryGraph.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/runGenerateTrajectoryGraph",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError occurred. otherwise r.json() will fail
        r.raise_for_status()
        # The index order relies on cells_id in an ascending form. The order is made in the R part.
        result = r.json()
        raise_if_error(result)

        return self._format_result(result)
