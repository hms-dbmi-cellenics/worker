import json
import os

import pytest
from exceptions import PythonWorkerException
from worker.helpers.get_diff_expr_cellsets import get_diff_expr_cellsets


class TestGetDiffExprCellSets:
    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "MockCellSet.json")) as f:
            self.cellsets = json.load(f)

    def test_should_throw_error_if_cell_sets_is_empty(self):
        basis_name = "patient-a"
        first_cell_set_name = "louvain-1"
        second_cell_set_name = "louvain-2"

        with pytest.raises(PythonWorkerException):
            get_diff_expr_cellsets(
                basis_name, first_cell_set_name, second_cell_set_name, self.cellsets
            )
