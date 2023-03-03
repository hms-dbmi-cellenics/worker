memoisedGetTopMarkerGenes <- memoise(
  getTopMarkerGenes,
  envir = .GlobalEnv,
  # cache_mem doesn't work because each request is run on a different process
  # so they don't share memory
  cache = cachem::cache_disk(
    dir="cache_marker_genes",
    max_size = 1024 * 1024^2,
    destroy_on_finalize = FALSE
  ),
  # Ignore scdata changing (its size makes it a bad idea to hash) use cleanup_cache instead
  # Ignore cell_sets, we use cell_sets_keys for caching
  omit_args = c("data", "cellSetsKeys")
)

hola <- function(a, b) {
  message("IMACTUALLYRUNNINGWTF")
  return("HOALHOLA")
}

memoisedHola <- memoise(
  hola,
  envir = .GlobalEnv,
  cache = cachem::cache_mem(max_size = 1024 * 1024^2)
  # Ignore scdata changing (its size makes it a bad idea to hash) use cleanup_cache instead
  # Ignore cell_sets, we use cell_sets_keys for caching
  # omit_args = c("cellSetsKeys")
)

# Cleans up all the caches that depend on the seurat object
# should be run whenever the seurat object changes
cleanup_cache <- function() {
  forget(memoisedGetTopMarkerGenes)
}

# create_cache() <- {
#   return(
#     memoisedGetTopMarkerGenes
#   );
# }