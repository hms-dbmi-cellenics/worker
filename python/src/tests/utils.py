def get_cell_ids(cell_class_key, cell_set_key, cell_sets):
    cell_class = next(cell_class for cell_class in cell_sets["cellSets"] if cell_class["key"] == cell_class_key)
    cell_ids = next(cell_set for cell_set in cell_class["children"] if cell_set["key"] == cell_set_key)["cellIds"]
    return cell_ids