import json
import os

import pytest
from exceptions import ErrorCodes, PythonWorkerException
from worker.helpers.get_diff_expr_cellsets import get_diff_expr_cellsets


class TestGetDiffExprCellSets:
    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "MockCellSet.json")) as f:
            self.cellsets = json.load(f)

        # @responses.activate

    # def test_throws_when_second_cellset_missing(self, mock_S3_get):
    #     MockS3Class.setResponse("one_set")
    #     with pytest.raises(Exception, match="fullfills the 2nd cell set"):
    #         DifferentialExpression(self.get_request())._format_request()

    # @responses.activate
    # def test_cells_in_sets_intersection_are_filtered_out(self, mock_S3_get):
    #     MockS3Class.setResponse("two_sets_intersected")

    #     request = DifferentialExpression(
    #         self.get_request(cellSet="cluster1", compareWith="cluster2")
    #     )._format_request()

    #     baseCells = request["baseCells"]
    #     backgroundCells = request["backgroundCells"]

    #     # Check 1 cell of each of the cell sets is left out
    #     assert len(baseCells) == len(backgroundCells) == 2

    #     # Check the cells that haven't been left out are
    #     # those that are not in the intersection of both sets
    #     assert len(set(baseCells).intersection(set(backgroundCells))) == 0

    def test_should_throw_error_if_1st_cell_sets_is_empty(self):
        basis_name = "patient-b"
        first_cell_set_name = "louvain-1"
        second_cell_set_name = "louvain-2"

        with pytest.raises(PythonWorkerException) as exc_info:
            get_diff_expr_cellsets(
                basis_name, first_cell_set_name, second_cell_set_name, self.cellsets
            )

        assert exc_info.value.args[0] == ErrorCodes.INVALID_INPUT
        assert exc_info.value.args[1] == "No cell id fullfills the 1st cell set."

    def test_should_throw_error_if_2nd_cell_sets_is_empty(self):
        basis_name = "patient-a"
        first_cell_set_name = "louvain-1"
        second_cell_set_name = "louvain-2"

        with pytest.raises(PythonWorkerException) as exc_info:
            get_diff_expr_cellsets(
                basis_name, first_cell_set_name, second_cell_set_name, self.cellsets
            )

        assert exc_info.value.args[0] == ErrorCodes.INVALID_INPUT
        assert exc_info.value.args[1] == "No cell id fullfills the 2nd cell set."
