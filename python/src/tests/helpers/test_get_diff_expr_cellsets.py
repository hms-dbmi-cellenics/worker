import json
import os

import pytest
from worker_status_codes import INVALID_INPUT
from exceptions import PythonWorkerException
from worker.helpers.get_diff_expr_cellsets import get_diff_expr_cellsets


class TestGetDiffExprCellSets:
    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "MockCellSet.json")) as f:
            cell_sets = json.load(f)
            self.cellsets = cell_sets["cellSets"]

    def test_should_throw_error_if_1st_cell_sets_is_empty(self):
        basis_name = "patient-b"
        first_cell_set_name = "louvain-1"
        second_cell_set_name = "louvain-2"

        with pytest.raises(PythonWorkerException) as exception_info:
            get_diff_expr_cellsets(
                basis_name, first_cell_set_name, second_cell_set_name, self.cellsets
            )

        assert exception_info.value.args[0] == INVALID_INPUT
        assert exception_info.value.args[1] == "No cell id fullfills the 1st cell set."

    def test_should_throw_error_if_2nd_cell_sets_is_empty(self):
        basis_name = "patient-a"
        first_cell_set_name = "louvain-1"
        second_cell_set_name = "louvain-2"

        with pytest.raises(PythonWorkerException) as exception_info:
            get_diff_expr_cellsets(
                basis_name, first_cell_set_name, second_cell_set_name, self.cellsets
            )

        assert exception_info.value.args[0] == INVALID_INPUT
        assert exception_info.value.args[1] == "No cell id fullfills the 2nd cell set."
