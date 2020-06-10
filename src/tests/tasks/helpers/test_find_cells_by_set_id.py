from tasks.helpers.find_cells_by_set_id import find_cells_by_set_id


class TestFindCellsBySetID:
    def test_empty_cell_set_returns_no_cells(self):
        assert find_cells_by_set_id("", []) == []

    def test_empty_cell_set_returns_appropriate_results(self):
        haystack = [
            {"cellIds": [1, 2, 3], "key": "asd"},
            {"cellIds": [4, 5, 6], "key": "fgh"},
        ]

        assert find_cells_by_set_id("asd", haystack) == [1, 2, 3]
        assert find_cells_by_set_id("fgh", haystack) == [4, 5, 6]

    def test_empty_cell_set_returns_appropriate_nested_results(self):
        haystack = [
            {"key": "asd", "children": [{"key": "fgh", "cellIds": [1, 2, 3]}]},
            {"cellIds": [4, 5, 6], "key": "ijk"},
        ]

        assert find_cells_by_set_id("fgh", haystack) == [1, 2, 3]
