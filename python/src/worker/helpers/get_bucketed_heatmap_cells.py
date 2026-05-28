from functools import reduce

def get_cell_class_ids(key, cell_sets):
  """Get all cell IDs from a cell class (root node)."""
  children = next(cell_class for cell_class in cell_sets if cell_class["key"] == key)["children"]
  
  cell_ids_set = set()
  for child in children:
    cell_ids_set.update(child["cellIds"])
  
  return cell_ids_set


def get_bucketed_heatmap_cells(selected_cell_set, grouped_tracks, cell_sets):
  """
  Downsample cells by grouping into buckets via cartesian product,
  then capping each bucket to 1000 cells.
  
  This is optimized for large datasets:
  - No hidden sets filtering (client handles that)
  - No selectedPoints filtering (client handles that)
  - Buckets capped at 1000 cells each (not proportional)
  
  Args:
    selected_cell_set: Key of the primary cell set (e.g., "louvain")
    grouped_tracks: Array of cell set keys for cartesian product
    cell_sets: Cell sets structure with hierarchy and properties
    
  Returns:
    List of cell IDs (up to 1000 per bucket)
  """
  # Get all cells from the selected cell set, intersected with louvain
  # to exclude filtered cell IDs (e.g. samples contains both filtered and unfiltered)
  enabled_cell_ids = get_cell_class_ids(selected_cell_set, cell_sets)
  filtered_cell_ids = get_cell_class_ids('louvain', cell_sets)
  enabled_cell_ids = enabled_cell_ids & filtered_cell_ids

  if not grouped_tracks or not enabled_cell_ids:
    return []
  
  # Build lookup for fast access
  cell_sets_by_key = {}
  for cell_class in cell_sets:
    cell_sets_by_key[cell_class["key"]] = cell_class
  
  def get_intersections(cell_ids, cell_class_key):
    """Get intersections of cell_ids with each child of the cell class."""
    children = cell_sets_by_key[cell_class_key]["children"]
    intersections = []
    
    for child in children:
      intersection = cell_ids & set(child["cellIds"])  # Fast set intersection
      if intersection:  # Only keep non-empty intersections
        intersections.append(intersection)
    
    return intersections
  
  def cartesian_product_intersection(buckets, cell_class_key):
    """Split each bucket by intersecting with children of cell_class_key."""
    new_buckets = []
    
    for bucket in buckets:
      intersections = get_intersections(bucket, cell_class_key)
      
      # Calculate leftover cells (cells not in any child)
      covered_cells = reduce(lambda acc, curr: acc | curr, intersections, set())
      leftover = bucket - covered_cells
      
      # Add all intersections
      new_buckets.extend(intersections)
      
      # Add leftover as its own bucket if non-empty
      if leftover:
        new_buckets.append(leftover)
    
    return new_buckets
  
  # Split into buckets by cartesian product of grouped tracks
  buckets = [enabled_cell_ids]
  
  for cell_class_key in grouped_tracks:
    buckets = cartesian_product_intersection(buckets, cell_class_key)
  
  # Cap each bucket to 1000 cells and combine
  import random
  result = []
  
  for bucket in buckets:
    bucket_list = list(bucket)
    sample_size = min(len(bucket_list), 1000)
    result.extend(random.sample(bucket_list, sample_size))
  
  return result