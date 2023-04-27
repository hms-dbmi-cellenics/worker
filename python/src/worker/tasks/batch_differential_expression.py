
import backoff
import requests
from aws_xray_sdk.core import xray_recorder
from exceptions import raise_if_error
import array
from ..tasks import Task
from ..result import Result
from ..config import config
from ..helpers.get_diff_expr_cellsets import get_diff_expr_cellsets
from ..helpers.s3 import get_cell_sets
import json
from exceptions import PythonWorkerException

class BatchDifferentialExpression(Task):
    def __init__(self, msg):
            super().__init__(msg)
            self.experiment_id = config.EXPERIMENT_ID

    def _format_result(self, results):
        # Return a list of formatted results.
        final_result=[]
        for data in results:
            final_result.append({"total": data["full_count"], "data": data["gene_results"]})
        return Result(final_result)
    
    def _format_request(self, base_cs, first_cs, second_cell_set_name, cell_sets):
        base_cells, background_cells = get_diff_expr_cellsets(
            str(base_cs), str(first_cs), second_cell_set_name, cell_sets
        )
        request = {
            "baseCells": [int(x) for x in base_cells],
            "backgroundCells": [int(x) for x in background_cells],
            "genesOnly": self.task_def.get("genesOnly", False),
            "comparisonType": self.task_def.get("comparisonType", "within"),
        }
        return request

    @xray_recorder.capture("DifferentialExpression.compute")
    @backoff.on_exception(
        backoff.expo, requests.exceptions.RequestException, max_time=30
    )

    def compute(self):
        # get cell sets from database
        cell_sets = get_cell_sets(self.experiment_id)
        first_cell_set_name = self.task_def["cellSet"]
        second_cell_set_name = self.task_def["compareWith"]
        basis = self.task_def["basis"]
        
        # either basis or first_cell_set are arrays, depending on what operation the user chose
        if len(basis) == 1:
            cell_sets_list = [(basis[0], cs) for cs in first_cell_set_name]
        else:
            cell_sets_list = [(b, first_cell_set_name[0]) for b in basis]

        responses_list = []

        for base_cs, first_cs in cell_sets_list:
            try:
                request = self._format_request(base_cs, first_cs, second_cell_set_name, cell_sets)
                 # send request to r worker
                response = requests.post(
                    f"{config.R_WORKER_URL}/v0/DifferentialExpression",
                    headers={"content-type": "application/json"},
                    data=json.dumps(request),
                )
                response.raise_for_status()
                result = response.json()
                raise_if_error(result)

                data = result.get("data")
            except PythonWorkerException:
                print(f"Cannot run Differential Expression for the current comparison, skipping...")
                data = {'full_count': 0, 'gene_results': 'No cell id fulfills the 1st cell set.'}

            responses_list.append(data)

        return self._format_result(responses_list)
