import json
import os

import pytest

from worker.helpers.find_cell_ids_in_same_hierarchy import (
    find_all_cell_ids_in_cell_sets,
    find_cell_ids_in_same_hierarchy,
)


class TestFindCellIdsInSameHierarchy:
    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "FindCell.json")) as f:
            self.cellsets = json.load(f)

    def test_empty_cell_set_returns_no_cells_hierchy(self):
        assert find_cell_ids_in_same_hierarchy("", []) == []

    def test_empty_cell_set_returns_no_cells(self):
        assert find_all_cell_ids_in_cell_sets([]) == []

    def test_empty_cell_set_returns_appropriate_results_hierarchy(self):
        assert find_cell_ids_in_same_hierarchy("louvain-11", self.cellsets) == [
            2,
            3,
            26,
            1094,
        ]
        assert find_cell_ids_in_same_hierarchy("condition-control", self.cellsets) == [
            4,
            18,
            1110,
        ]

    def test_empty_cell_set_returns_appropriate_results(self):
        assert find_all_cell_ids_in_cell_sets(self.cellsets) == [
            2,
            3,
            26,
            1094,
            86,
            1105,
            0,
            21,
            29,
            31,
            1113,
            4,
            18,
            1110,
        ]
