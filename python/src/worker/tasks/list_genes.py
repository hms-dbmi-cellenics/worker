import json

import backoff
import pandas as pd
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.remove_regex import remove_regex
from ..result import Result
from ..tasks import Task


class ListGenes(Task):
    def _format_result(self, result, total):
        # convert result to list of row dicts
        result = result.to_dict(orient="records")

        # Return a list of formatted results.
        return Result({"total": total, "rows": result})

    @xray_recorder.capture("ListGenes.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
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

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/listGenes",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise
        #  r.json() will fail
        r.raise_for_status()
        r = r.json()
        total = r["full_count"]
        result = pd.DataFrame.from_dict(r["gene_results"])
        return self._format_result(result, total=total)
