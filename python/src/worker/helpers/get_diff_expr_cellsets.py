from exceptions import PythonWorkerException
from worker_status_codes import INVALID_INPUT

from .find_cell_ids_in_same_hierarchy import (find_all_cell_ids_in_cell_sets,
                                              find_cell_ids_in_same_hierarchy)
from .find_cells_by_set_id import find_cells_by_set_id


def get_diff_expr_cellsets(
    basis_name,
    first_cell_set_name,
    second_cell_set_name,
    all_cell_sets,
):
    # Check if the comparsion is between all the cells or within a cluster
    if not basis_name or ("all" in basis_name.lower()):
        # In not filtering by a cluster we leave the set empty
        filtered_set = set()
    else:
        # if filtering by a cluster we keep the cell ids in a set
        filtered_set = set(find_cells_by_set_id(basis_name, all_cell_sets))

    # mark cells of first set
    first_cell_set = set(find_cells_by_set_id(first_cell_set_name, all_cell_sets))

    # mark cells of second set
    # check if the second set is composed by the "All other cells"
    if second_cell_set_name == "background" or "all" in second_cell_set_name.lower():
        # Retrieve all cells (not necessarily at the same hierarchy level)
        complete_cell_set = set(find_all_cell_ids_in_cell_sets(all_cell_sets))
        # Filter with those that are not in the first cell set
        second_cell_set = complete_cell_set.difference(first_cell_set)
    else:
        # In the case that we compare with specific cell set, we just look for
        #  the cell directly
        second_cell_set = set(
            get_cells_in_set(first_cell_set_name, second_cell_set_name, all_cell_sets)
        )
        # Check any possible intersect cells
        inter_cell_set = first_cell_set.intersection(second_cell_set)

        first_cell_set = first_cell_set.difference(inter_cell_set)
        second_cell_set = second_cell_set.difference(inter_cell_set)

    # Keep only cells that are on the filtered basis (if not in "All" analysis)
    if len(filtered_set) > 0:
        second_cell_set = second_cell_set.intersection(filtered_set)
        first_cell_set = first_cell_set.intersection(filtered_set)

    # Check if the first cell set is empty
    if len(first_cell_set) == 0:
        raise PythonWorkerException(
            INVALID_INPUT, "No cell id fullfills the 1st cell set."
        )

    # Check if the second cell set is empty
    if len(second_cell_set) == 0:
        raise PythonWorkerException(
            INVALID_INPUT, "No cell id fullfills the 2nd cell set."
        )

    return first_cell_set, second_cell_set


# Get cells for the cell set.
def get_cells_in_set(first_cell_set_name, second_cell_set_name, all_cell_sets):
    cells = []

    # If "rest", then get all cells in the same hierarchy as the first cell set
    #  that arent part of "first"
    if "rest" in second_cell_set_name.lower():
        cells = find_cell_ids_in_same_hierarchy(first_cell_set_name, all_cell_sets)
    else:
        cells = find_cells_by_set_id(second_cell_set_name, all_cell_sets)

    return cells
