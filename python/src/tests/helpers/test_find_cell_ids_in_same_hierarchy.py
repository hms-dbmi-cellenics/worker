import json

from worker.helpers.find_cell_ids_in_same_hierarchy import (
    find_all_cell_ids_in_cell_sets,
    find_cell_ids_in_same_hierarchy,
)

cells_data_file = "tests/data/FindCell.json"


class TestFindCellIdsInSameHierarchy:
    def test_empty_cell_set_returns_no_cells_hierchy(self):
        assert find_cell_ids_in_same_hierarchy("", []) == []

    def test_empty_cell_set_returns_no_cells(self):
        assert find_all_cell_ids_in_cell_sets([]) == []

    def test_empty_cell_set_returns_appropriate_results_hierarchy(self):
        with open(cells_data_file) as f:
            haystack = json.load(f)
        print(find_cell_ids_in_same_hierarchy("condition-control", haystack))
        assert find_cell_ids_in_same_hierarchy("louvain-11", haystack) == [
            2,
            3,
            26,
            1094,
        ]
        assert find_cell_ids_in_same_hierarchy("condition-control", haystack) == [
            4,
            18,
            1110,
        ]

    def test_empty_cell_set_returns_appropriate_results(self):
        with open(cells_data_file) as f:
            haystack = json.load(f)
        print(find_all_cell_ids_in_cell_sets(haystack))
        assert find_all_cell_ids_in_cell_sets(haystack) == [
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
