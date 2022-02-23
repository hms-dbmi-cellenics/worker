import json
import os

import pytest
from worker.helpers.find_cell_ids_in_same_hierarchy import (
    find_all_cell_ids_in_cell_sets, find_cell_ids_in_same_hierarchy)


class TestFindCellIdsInSameHierarchy:
    @pytest.fixture(autouse=True)
    def load_cellsets(self):
        with open(os.path.join("tests/data", "MockCellSet.json")) as f:
            self.cellsets = json.load(f)

    def test_empty_cell_set_returns_no_cells_hierchy(self):
        assert find_cell_ids_in_same_hierarchy("", []) == []

    def test_empty_cell_set_returns_no_cells(self):
        assert find_all_cell_ids_in_cell_sets([]) == []

    def test_empty_cell_set_returns_appropriate_results_hierarchy(self):
        print(find_cell_ids_in_same_hierarchy("condition-control", self.cellsets))
        assert find_cell_ids_in_same_hierarchy("louvain-2", self.cellsets) == [
            1,
            2,
            3,
            4,
            5,
        ]
        assert find_cell_ids_in_same_hierarchy("condition-control", self.cellsets) == [
            4,
            5,
            6,
        ]

    def test_empty_cell_set_returns_appropriate_results(self):
        print(find_all_cell_ids_in_cell_sets(self.cellsets))
        assert find_all_cell_ids_in_cell_sets(self.cellsets) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
