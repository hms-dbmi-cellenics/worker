def find_cells_by_set_id(needle, haystack):
    for cell_set in haystack:
        if cell_set["key"] == needle:
            return cell_set["cellIds"]

        if cell_set.get("children", None):
            result = find_cells_by_set_id(needle, cell_set["children"])

            if result:
                return result

    return []
