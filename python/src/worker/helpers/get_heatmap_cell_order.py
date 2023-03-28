from functools import reduce
import math
import random

from .cell_sets_dict import get_cell_sets_dict

def get_cell_class_ids(key, cell_sets):
  children = next(cell_class for cell_class in cell_sets if cell_class["key"] == key)["children"]
  cell_ids = map(lambda cell_set: cell_set["cellIds"], children)

  cell_ids_set = set()
  for i in cell_ids:
    cell_ids_set.update(i)

  return cell_ids_set

def get_heatmap_cell_order(selected_cell_set, grouped_tracks, selected_points, hidden_cell_set_keys, max_cells, cell_sets):
  cell_sets_by_key = get_cell_sets_dict(cell_sets)

  filtered_cell_ids = get_cell_class_ids('louvain', cell_sets)

  def get_cells(key, is_root_node=False):
    unfiltered_cell_ids = None

    if (is_root_node): 
      unfiltered_cell_ids = get_cell_class_ids(key, cell_sets)
    else: 
      unfiltered_cell_ids = cell_sets_by_key[key]["cellIds"]

    return filtered_cell_ids.intersection(set(unfiltered_cell_ids))

  def get_all_enabled_cell_ids():
    cell_ids = get_cells(selected_cell_set, is_root_node=True)

    if (selected_points != "All"):
      cell_set_key = selected_points.split('/')[1]
      cell_ids = cell_ids.intersection(get_cells(cell_set_key))
      
    for hidden_cell_set in hidden_cell_set_keys:
      cell_ids = cell_ids.difference(get_cells(hidden_cell_set))

    return cell_ids

  # Returns a list of sets, one for each cell set in the cell_class 
  # each of them is an intersection of each cell set with cell_ids
  def get_intersections(cell_ids, cell_class):
    children = next(curr_class for curr_class in cell_sets if curr_class["key"] == cell_class)["children"]
    cell_ids_by_set = map(lambda cell_set: cell_set["cellIds"], children)

    intersections = []
    for current_set_ids in cell_ids_by_set:

      current_intersection = set(current_set_ids).intersection(cell_ids)
      
      if (len(current_intersection) > 0):
        intersections.append(current_intersection)

    return intersections

  def cartesian_product_intersection(buckets, cell_class):
    new_buckets = []

    for bucket in buckets:
      intersections = get_intersections(bucket, cell_class)

      # The cells that werent part of any intersection are also added at the end
      leftover_cells = reduce(lambda acum, current: acum.difference(current), intersections, bucket)
      intersections.append(leftover_cells)
      
      for intersection in intersections:
        new_buckets.append(intersection)

    return new_buckets

  def split_by_cartesian_intersections(enabled_cell_ids):
    buckets = None
    size = None

    buckets = [enabled_cell_ids]

    for cell_class_key in grouped_tracks:
      buckets = cartesian_product_intersection(buckets, cell_class_key)

    # We need to calculate size at the end because we may have repeated cells
    # (due to group bys having the same cell in different groups)
    size = reduce(lambda acum, bucket: acum + len(bucket), buckets, 0)

    return buckets, size
  
  def downsample(buckets, amount_of_cells):
    downsampled_cell_ids = []

    # If we collected less than max_cells, then no need to downsample
    final_sample_size = min(amount_of_cells, max_cells)

    for bucket in buckets:
      sample_size = math.floor((len(bucket) / amount_of_cells) * final_sample_size)

      # Pick sample_size elements randomly
      sample = random.sample(list(bucket), sample_size)

      downsampled_cell_ids.extend(sample)

    return downsampled_cell_ids
  
  
  enabled_cell_ids = get_all_enabled_cell_ids()
  
  if (len(grouped_tracks) == 0 or len(enabled_cell_ids) == 0): 
    return []

  buckets, size = split_by_cartesian_intersections(enabled_cell_ids)

  return downsample(buckets, size)
