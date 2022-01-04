#
# getClusters
# returns the clusters in the shape of a dataframe with a clusters column,
# cell ids column and cell barcode as rownames.
#
# req$body has:
# type: can be "louvain"/"leiden"
# config:{
#          resolution: integer, range: 0 - 2
#         }
#
#
# We currently CANT support leiden, we need to discuss this in bioinformatics,
#  the algorithm is not working.
#
#' @export
#'
runClusters <- function(req, data) {
  resol <- req$body$config$resolution
  type <- req$body$type

  data <- getClusters(type, resol, data)
  res_col <- paste0(data@active.assay, "_snn_res.", toString(resol))
  # In the meta data slot the clustering is stored with the resolution
  # used to calculate it
  # RNA_snn_res.#resolution
  df <- data.frame(
    cluster = data@meta.data[, res_col],
    cell_ids = data@meta.data$cells_id
  )
  # get the cell barcodes as rownames
  rownames(df) <- rownames(data@meta.data)
  return(df)
}
