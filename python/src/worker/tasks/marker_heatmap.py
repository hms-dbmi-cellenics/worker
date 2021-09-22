import json

import backoff
import numpy as np
import requests
from aws_xray_sdk.core import xray_recorder

from ..config import config
from ..helpers.s3 import get_cell_sets
from ..result import Result
from ..tasks import Task


class MarkerHeatmap(Task):
    def __init__(self, msg):
        super().__init__(msg)
        self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, result):
        # JSONify result.
        result = json.dumps(result)
        # Return a list of formatted results.
        return [Result(result)]

    @xray_recorder.capture("GeneExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )
    def compute(self):
        request = {"nGenes": self.task_def["nGenes"]}
        
        cellSetKey = self.task_def["cellSetKey"]

        cellSets = get_cell_sets(self.experiment_id)

        for set in cellSets:
            if(set["key"]==cellSetKey):
                cellSets = set
                break
        
        request["cellSets"] = cellSets

        r = requests.post(
            f"{config.R_WORKER_URL}/v0/runMarkerHeatmap",
            headers={"content-type": "application/json"},
            data=json.dumps(request),
        )
        # raise an exception if an HTTPError if one occurred because otherwise r.json() will fail
        r.raise_for_status()
        resultR = r.json()
        truncatedR = resultR["truncatedExpression"]
        resultR = resultR["rawExpression"]
        result = {}
        data = {}
        order = []
        if not len(resultR):
            result[genes[0]] = {
                "error": 404,
                "message": "Gene {} not found!".format(genes[0]),
            }
        else:
            for gene in resultR.keys():
                view = resultR[gene]
                # can't do summary stats on list with None's
                # casting to np array replaces None with np.nan
                viewnp = np.array(view, dtype=np.float)
                # This is not necessary and is also costly, but I leave it commented as a reminder
                # that this object has integer zeros and floating point for n!=0.
                # expression = [float(item) for item in view]
                mean = float(np.nanmean(viewnp))
                stdev = float(np.nanstd(viewnp))
                data[gene] = {"truncatedExpression": {}, "rawExpression": {}}
                data[gene]["rawExpression"] = {
                    "mean": mean,
                    "stdev": stdev,
                    "expression": view,
                }

                viewTr = truncatedR[gene]
                viewnpTr = np.array(viewTr, dtype=np.float)
                minimum = float(np.nanmin(viewnpTr))
                maximum = float(np.nanmax(viewnpTr))
                data[gene]["truncatedExpression"] = {
                                    "min": minimum,
                                    "max": maximum,
                                    "expression": viewTr,
                                }    
                order.append(gene)  
        result["data"] = data
        result["order"] = order   
        return self._format_result(result)
