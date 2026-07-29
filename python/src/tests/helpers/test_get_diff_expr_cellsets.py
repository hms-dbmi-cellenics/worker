import copy
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

    # Cell set keys can be derived from user facing names, so they can contain
    # "all" or "rest" as a substring (e.g. a cell type called "Small
    # intestine", a group called "Forest")
    def add_cellset(self, cell_class_key, key, cell_ids):
        cellsets = copy.deepcopy(self.cellsets)
        cell_class = next(cs for cs in cellsets if cs["key"] == cell_class_key)
        cell_class["children"].append(
            {"key": key, "name": key, "color": "#000000", "cellIds": cell_ids}
        )

        return cellsets

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

    def test_should_filter_by_a_basis_that_contains_all_in_its_name(self):
        # "Small" contains "all", but this is a cell set like any other
        cellsets = self.add_cellset("louvain", "louvain-Small intestine", [1, 2, 4])

        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            "louvain-Small intestine",
            "condition-control",
            "condition-treated",
            cellsets,
        )

        assert first_cell_set == {1, 2}
        assert second_cell_set == {4}

    def test_should_compare_with_a_cell_set_that_contains_all_in_its_name(self):
        cellsets = self.add_cellset("condition", "condition-Small intestine", [7, 8])

        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            None, "condition-control", "condition-Small intestine", cellsets
        )

        assert first_cell_set == {1, 2, 3}
        assert second_cell_set == {7, 8}

    def test_should_compare_with_a_cell_set_that_contains_rest_in_its_name(self):
        cellsets = self.add_cellset("condition", "condition-Forest", [7, 8])

        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            None, "condition-control", "condition-Forest", cellsets
        )

        assert first_cell_set == {1, 2, 3}
        assert second_cell_set == {7, 8}

    def test_should_compare_with_the_rest_of_the_cell_sets(self):
        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            None, "condition-control", "rest", self.cellsets
        )

        assert first_cell_set == {1, 2, 3}
        assert second_cell_set == {4, 5, 6}

    def test_should_compare_with_all_the_other_cells(self):
        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            None, "condition-control", "background", self.cellsets
        )

        assert first_cell_set == {1, 2, 3}
        assert second_cell_set == {4, 5, 6, 7, 8, 9, 10}

    def test_should_not_filter_by_basis_when_comparing_all_the_cells(self):
        first_cell_set, second_cell_set = get_diff_expr_cellsets(
            "all", "condition-control", "condition-treated", self.cellsets
        )

        assert first_cell_set == {1, 2, 3}
        assert second_cell_set == {4, 5, 6}
