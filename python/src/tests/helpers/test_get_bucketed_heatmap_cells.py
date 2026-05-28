import pytest

from worker.helpers.get_bucketed_heatmap_cells import (
    get_bucketed_heatmap_cells,
    get_cell_class_ids,
)

CELL_SETS = [
    {
        "key": "louvain",
        "children": [
            {"key": "louvain-0", "cellIds": [0, 1, 2]},
            {"key": "louvain-1", "cellIds": [3, 4, 5]},
        ],
    },
    {
        "key": "sample",
        "children": [
            {"key": "sample-a", "cellIds": [0, 1, 3, 4]},
            {"key": "sample-b", "cellIds": [2, 5]},
        ],
    },
    {
        "key": "scratchpad",
        "children": [],
    },
]


class TestGetCellClassIds:
    def test_returns_ids_from_all_children(self):
        ids = get_cell_class_ids("louvain", CELL_SETS)
        assert ids == {0, 1, 2, 3, 4, 5}

    def test_returns_ids_from_single_child(self):
        ids = get_cell_class_ids("sample", CELL_SETS)
        assert ids == {0, 1, 2, 3, 4, 5}

    def test_empty_children_returns_empty_set(self):
        ids = get_cell_class_ids("scratchpad", CELL_SETS)
        assert ids == set()


class TestGetBucketedHeatmapCells:
    def test_no_grouped_tracks_returns_empty(self):
        result = get_bucketed_heatmap_cells("louvain", [], CELL_SETS)
        assert result == []

    def test_empty_selected_cell_set_returns_empty(self):
        empty_cell_sets = [
            {"key": "louvain", "children": []},
            {"key": "sample", "children": []},
        ]
        result = get_bucketed_heatmap_cells(
            "louvain", ["sample"], empty_cell_sets
        )
        assert result == []

    def test_returns_subset_of_enabled_cells(self):
        result = get_bucketed_heatmap_cells(
            "louvain", ["sample"], CELL_SETS
        )
        # All returned IDs must come from louvain (the enabled set)
        louvain_ids = {0, 1, 2, 3, 4, 5}
        assert all(cell_id in louvain_ids for cell_id in result)

    def test_result_length_at_most_1000_per_bucket(self):
        # Build a large cell set with 2000 cells in one bucket
        large_cell_sets = [
            {
                "key": "louvain",
                "children": [{"key": "louvain-0", "cellIds": list(range(2000))}],
            },
            {
                "key": "sample",
                "children": [{"key": "sample-a", "cellIds": list(range(2000))}],
            },
        ]
        result = get_bucketed_heatmap_cells(
            "louvain", ["sample"], large_cell_sets
        )
        assert len(result) <= 1000

    def test_single_bucket_contains_correct_cells(self):
        # One track, two samples → two buckets (intersect louvain with samples)
        result = get_bucketed_heatmap_cells(
            "louvain", ["sample"], CELL_SETS
        )
        result_set = set(result)
        # IDs in result must be a subset of louvain IDs
        assert result_set <= {0, 1, 2, 3, 4, 5}

    def test_multiple_grouped_tracks_produce_cartesian_buckets(self):
        cell_sets = [
            {
                "key": "louvain",
                "children": [
                    {"key": "louvain-0", "cellIds": [0, 1, 2, 3]},
                ],
            },
            {
                "key": "sample",
                "children": [
                    {"key": "sample-a", "cellIds": [0, 1]},
                    {"key": "sample-b", "cellIds": [2, 3]},
                ],
            },
            {
                "key": "condition",
                "children": [
                    {"key": "cond-ctrl", "cellIds": [0, 2]},
                    {"key": "cond-stim", "cellIds": [1, 3]},
                ],
            },
        ]
        result = get_bucketed_heatmap_cells(
            "louvain", ["sample", "condition"], cell_sets
        )
        result_set = set(result)
        # All cells should still be from louvain
        assert result_set <= {0, 1, 2, 3}
        # All original cells should be represented (small enough for 1000 cap)
        assert result_set == {0, 1, 2, 3}
