# Convert the cell sets object into a dictionary
# e.g. { 'louvain' : ..., 'louvain-1', ... }
# cell_sets is a cell sets object
def get_cell_sets_dict(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        if not cell_sets_dict.get(cell_class["key"]):
          cell_sets_dict[cell_class["key"]] = {
              "key": cell_class["key"],
              "name": cell_class["name"],
              "rootNode": True,
              "childrenKeys": [],
          }

        for cell_set in cell_class["children"]:
            cell_sets_dict[cell_class["key"]]["childrenKeys"].append(cell_set["key"])
            cell_sets_dict[cell_set["key"]] = cell_set
            cell_sets_dict[cell_set["key"]]["rootNode"] = False

    return cell_sets_dict

# Get all cell sets that match the subset_keys
# subset_keys: Array of cell set keys, e.g. ['louvain', 'louvain-1', ...]
# cell_sets_dict : Dictionary of cell sets. output of get_cell_sets_dict()
def subset_cell_sets_dict(subset_keys, cell_sets_dict):
    cell_ids = set()

    for subset_key in subset_keys:

        # If cell class, get cell ids from child keys
        if cell_sets_dict[subset_key]["rootNode"]:
            cell_ids = cell_ids.union(
                subset_cell_sets_dict(cell_sets_dict[subset_key]["childrenKeys"], cell_sets_dict)
            )
        else:
          cell_ids = cell_ids.union(set(cell_sets_dict[subset_key]['cellIds']))

    return cell_ids


# Convert the cell sets object into a dictionary
# that is easily parsable in R
def get_cell_sets_dict_for_r(cell_sets):
    cell_sets_dict = {}

    for cell_class in cell_sets:
        cell_sets_dict[cell_class["key"]] = cell_class

    return cell_sets_dict