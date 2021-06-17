import backoff
import pandas as pd
from result import Result
import requests
from config import get_config
from tasks import Task
import json
from aws_xray_sdk.core import xray_recorder

config = get_config()


class ListGenes(Task):
    def __init__(self, msg):
        self.task_def = msg["body"]

    def _format_result(self, result, total):
        # convert result to list of row dicts
        result = result.to_dict(orient="records")

        # JSONify result.
        result = json.dumps({"total": total, "rows": result})

        # Return a list of formatted results.
        return [Result(result)]

    @xray_recorder.capture('ListGenes.compute')
    @backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=30)
    def compute(self):
        request = self.task_def
        #
        # Remove potentialy hazardous characters from the names filter
        # Leaving this'd leave us open to a DDOS attack in the form of a very time consuming regex.
        # More info here: https://owasp.org/www-community/attacks/Regular_expression_Denial_of_Service_-_ReDoS
        #
        # The symbols ^ and $ can't be removed because they are the way the UI indicates the type of search the user is performing.
        #
        if "geneNamesFilter" in request:
            geneFilter = request["geneNamesFilter"]
            regexChars = "{}|()?¿*+|\/.<>"
            for char in regexChars:
                geneFilter = geneFilter.replace(char, "")
            request["geneNamesFilter"] = geneFilter
        r = requests.post(
            f"{config.R_WORKER_URL}/v0/listGenes",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )

        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        resR = r.json()

        # Convert to dataframe to prepare the data for the UI
        resR = pd.DataFrame(resR)
        total = 0
        if len(resR) > 0:
            total = resR["full_count"][0]
        resR = resR.drop("full_count", axis=1)
        # total returns numpy int64, convert to integer
        # for serialization to JSON
        total = int(total)
        return self._format_result(resR, total=total)