import json

import backoff
import pandas as pd
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error

from ..config import config
from ..helpers.remove_regex import remove_regex
from ..result import Result
from ..tasks import Task

class ListGenes(Task):
    def _format_result(self, result, total):
        # Return a list of formatted results.
        return Result({"total": total, "rows": result})

    def _format_request(self):
        request = self.task_def
        #
        # Remove potentialy hazardous characters from the names filter
        # Leaving this'd leave us open to a DDOS attack in the form of a very time
        #  consuming regex.
        # More info here:
        # https://owasp.org/www-community/attacks/Regular_expression_Denial_of_Service_-_ReDoS
        #
        # The symbols ^ and $ can't be removed because they are the way the UI
        # indicates the type of search the user is performing.
        #
        if "geneNamesFilter" in request:
            gene_filter = request["geneNamesFilter"]
            request["geneNamesFilter"] = remove_regex(gene_filter)

        return request

    @xray_recorder.capture("ListGenes.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = self._format_request()

        response = requests.post(
            f"{config.R_WORKER_URL}/v0/listGenes",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        response.raise_for_status()
        result = response.json()
        raise_if_error(result)

        data = result.get("data")

        total = data["full_count"]
        result = data["gene_results"]

        return self._format_result(result, total)
