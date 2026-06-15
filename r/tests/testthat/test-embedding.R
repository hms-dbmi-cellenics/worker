library("mockery")

mock_req <- function() {
  req <- list(
    body = list(
      type = "umap",
      config = list(minimumDistance = 0.1, distanceMetric = "cosine"),
      use_saved = FALSE
    )
  )
}

test_that("TSNE embedding works", {

  mock_RunTSNE <- function(config, method, reduction_type, num_pcs, data) {
    Seurat::RunTSNE(
      data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate,
      check_duplicates = FALSE
    )
  }

  reduction_method = 'tsne'
  stub(runEmbedding, 'getEmbedding', mock_RunTSNE)

  data <- suppressWarnings(mock_scdata())
  req <- list(
    body = list(
      type = reduction_method,
      config = list(perplexity = 10, learningRate = 100),
      use_saved = FALSE
    )
  )

  res <- runEmbedding(req, data)
  expected_res <- as.data.frame(Seurat::Embeddings(data)[,1:2])

  # Expect all cells to be in embedding
  expect_equal(length(res), length(expected_res$PC_1))
})

test_that("UMAP embedding works", {

  data <- mock_scdata(with_umap = TRUE)
  req <- list(
    body = list(
      type = "umap",
      config = list(minimumDistance = 0.1, distanceMetric = "cosine"),
      use_saved = FALSE
    )
  )

  res <- runEmbedding(req, data)
  expected_res <- as.data.frame(Seurat::Embeddings(data)[, 1:2])

  # Expect all cells to be in embedding
  expect_equal(length(res), length(expected_res$PC_1))
})

test_that("UMAP embedding works with bpcells", {

  data <- mock_scdata(use_bpcells = TRUE)
  req <- list(
    body = list(
      type = "umap",
      config = list(minimumDistance = 0.1, distanceMetric = "cosine"),
      use_saved = FALSE
    )
  )

  expect_no_error(runEmbedding(req, data))
})

test_that("Projection onto sketched embedding works", {

  data <- mock_scdata(use_bpcells = TRUE, nreps = 10)
  data <- suppressWarnings(mock_sketch(data))
  req <- list(
    body = list(
      type = "umap",
      config = list(minimumDistance = 0.1, distanceMetric = "cosine"),
      use_saved = FALSE
    )
  )

  expect_no_error(runEmbedding(req, data))
})

test_that("RunTSNE uses the correct params", {

  mock_RunTSNE <- mock(TRUE)

  stub(getEmbedding, "Seurat::RunTSNE", mock_RunTSNE)

  data <- suppressWarnings(mock_scdata())
  config <- list(perplexity = 10, learningRate = 100)
  reduction_type <- "pca"
  method <- "tsne"
  num_pcs <- 1

  res <- getEmbedding(config, method, reduction_type, num_pcs, data)

  # Check that tsne is called using the correct parameters
  expect_equal(length(mock_RunTSNE), 1)
  args <- mock_args(mock_RunTSNE)

  expect_equal(
    args[[1]],
    list(
      data,
      reduction = reduction_type,
      dims = 1:num_pcs,
      perplexity = config$perplexity,
      learning.rate = config$learningRate
    )
  )

})

test_that("assignEmbedding assigns embedding correctly for UMAP", {

  # given an embedding, which is ordered by cell id
  num_cells = 80
  mock_embedding <- list()

  for(i in 1:(num_cells * 2)) {
    if(i %% 2 != 0) {
      mock_embedding[[i]] <- c(i + 2, i + 2)
    } else {
      mock_embedding[[i]] <- c(NA, NA)
    }
  }

  # and a Seurat object with shuffled cell_ids order
  mock_seurat_object <- mock_scdata()
  old_embedding <- Seurat::Embeddings(mock_seurat_object)

  # cells_id 0 corresponds to embedding[[1]]
  mock_seurat_object$cells_id <- seq((num_cells * 2) - 2, 0, -2)

  # assigning embedding
  mock_seurat_object <- assignEmbedding(mock_embedding, mock_seurat_object)
  new_embedding <- Seurat::Embeddings(mock_seurat_object, reduction = "umap")

  # check that assignEmbedding replaces the embedding
  expect_false(identical(old_embedding, new_embedding))

  # check that assignEmbedding correctly orders the embedding
  expect_equal(rownames(new_embedding), colnames(mock_seurat_object))

  # check that assignEmbedding replaces with the right value
  test_cell_id <- 52

  barcode <- rownames(mock_seurat_object@meta.data[which(mock_seurat_object@meta.data$cells_id == test_cell_id), ])
  resulting_embedding <- new_embedding[barcode, ]

  # Add 1 to test_cell_id because cell_id 0 corresponds to embedding [[1]]
  expect_equal(
    unname(resulting_embedding),
    mock_embedding[[test_cell_id + 1]]
  )
})


test_that("assignEmbedding assigns embedding correctly for tSNE", {

  # given an embedding, which is ordered by cell id
  num_cells = 80
  mock_embedding <- list()

  for(i in 1:(num_cells * 2)) {
    if(i %% 2 != 0) {
      mock_embedding[[i]] <- c(i + 2, i + 2)
    } else {
      mock_embedding[[i]] <- c(NA, NA)
    }
  }

  # and a Seurat object with shuffled cell_ids order
  mock_seurat_object <- mock_scdata()
  old_embedding <- Seurat::Embeddings(mock_seurat_object)

  # cells_id 0 corresponds to embedding[[1]]
  mock_seurat_object$cells_id <- seq((num_cells * 2) - 2, 0, -2)

  # assigning embedding
  mock_seurat_object <- assignEmbedding(mock_embedding, mock_seurat_object, reduction_method = "tsne")
  new_embedding <- Seurat::Embeddings(mock_seurat_object, reduction = "tsne")

  # check that assignEmbedding replaces the embedding
  expect_false(identical(old_embedding, new_embedding))

  # check that assignEmbedding correctly orders the embedding
  expect_equal(rownames(new_embedding), colnames(mock_seurat_object))

  # check that assignEmbedding replaces with the right value
  test_cell_id <- 52

  barcode <- rownames(mock_seurat_object@meta.data[which(mock_seurat_object@meta.data$cells_id == test_cell_id), ])
  resulting_embedding <- new_embedding[barcode, ]

  # Add 1 to test_cell_id because cell_id 0 corresponds to embedding [[1]]
  expect_equal(unname(resulting_embedding), mock_embedding[[test_cell_id + 1]])
})

test_that("can request saved embedding result", {
  data <- suppressWarnings(mock_scdata())
  req <- mock_req()
  req$body$type <- "pca"
  req$body$use_saved <- TRUE

  res <- runEmbedding(req, data)

  expected_res <- as.data.frame(Seurat::Embeddings(data)[,1:2])

  expected_res <- expected_res |>
    as.data.frame() |>
    dplyr::rowwise() |>
    dplyr::mutate(PCS = list(c(PC_1, PC_2)))

  expect_equal(res,expected_res$PCS)
})

# ---------------------------------------------------------------------------
# Spatial embedding helpers (init-spatial)
# ---------------------------------------------------------------------------

# Minimal S4 mocks that mimic the slots/accessors used by the spatial helpers
# in embedding.R, so we don't need to construct a full VisiumV2/FOV object.
setClass("MockSeg", representation(sf.data = "data.frame"))
setClass(
  "MockImg",
  representation(image = "ANY", seg = "MockSeg", sf = "list")
)
# img[["segmentations"]] returns the segmentation object
setMethod("[[", "MockImg", function(x, i, ...) x@seg)

mock_img <- function(dims = c(500, 500, 3),
                     sf.data = data.frame(
                       cell = c("c1", "c2"),
                       x = c(2, 4),
                       y = c(3, 5)
                     ),
                     scalefactors = list(lowres = 0.1, hires = 0.5)) {
  new(
    "MockImg",
    image = array(0, dim = dims),
    seg = new("MockSeg", sf.data = sf.data),
    sf = scalefactors
  )
}

test_that("get_img_scale returns lowres when an image dim is <= 600", {
  scdata <- list(fov1 = mock_img(dims = c(500, 500, 3)))
  expect_equal(get_img_scale("fov1", scdata), "lowres")
})

test_that("get_img_scale returns hires when both image dims are > 600", {
  scdata <- list(fov1 = mock_img(dims = c(2000, 2000, 3)))
  expect_equal(get_img_scale("fov1", scdata), "hires")
})

test_that("get_img_scale only considers the first two (spatial) dims", {
  # 3rd dim (channels) is small but must be ignored
  scdata <- list(fov1 = mock_img(dims = c(1000, 1000, 3)))
  expect_equal(get_img_scale("fov1", scdata), "hires")
})

test_that("get_tissue_coords uses the resolved scale and tags cells", {
  scdata <- list(fov1 = mock_img(dims = c(500, 500, 3)))

  fake_coords <- data.frame(
    x = c(1, 2),
    y = c(3, 4),
    row.names = c("c1", "c2")
  )
  mock_gtc <- mock(fake_coords)
  stub(get_tissue_coords, "SeuratObject::GetTissueCoordinates", mock_gtc)

  res <- get_tissue_coords("fov1", scdata)

  # cell column is added from rownames
  expect_equal(res$cell, c("c1", "c2"))
  # coordinates are passed through unrotated (no [, 2:1] flip)
  expect_equal(res$x, c(1, 2))
  expect_equal(res$y, c(3, 4))

  # GetTissueCoordinates called with the scale picked by get_img_scale (lowres)
  args <- mock_args(mock_gtc)[[1]]
  expect_equal(args$scale, "lowres")
})

test_that("get_polygon_coords scales x/y by the resolved scale factor", {
  img <- mock_img(
    dims = c(500, 500, 3),
    sf.data = data.frame(
      cell = c("c1", "c2"),
      x = c(2, 4),
      y = c(3, 5)
    ),
    scalefactors = list(lowres = 0.1, hires = 0.5)
  )
  scdata <- list(fov1 = img)

  # ScaleFactors(img) -> the named list of scale factors
  stub(get_polygon_coords, "Seurat::ScaleFactors", function(...) img@sf)

  res <- get_polygon_coords("fov1", scdata)

  # dims <= 600 -> lowres factor 0.1 applied to x and y
  expect_equal(res$x, c(2, 4) * 0.1)
  expect_equal(res$y, c(3, 5) * 0.1)
  # non-coordinate columns are untouched
  expect_equal(res$cell, c("c1", "c2"))
})

test_that("get_polygon_coords uses hires factor for large images", {
  img <- mock_img(
    dims = c(2000, 2000, 3),
    sf.data = data.frame(
      cell = c("c1", "c2"),
      x = c(2, 4),
      y = c(3, 5)
    ),
    scalefactors = list(lowres = 0.1, hires = 0.5)
  )
  scdata <- list(fov1 = img)

  stub(get_polygon_coords, "Seurat::ScaleFactors", function(...) img@sf)

  res <- get_polygon_coords("fov1", scdata)

  expect_equal(res$x, c(2, 4) * 0.5)
  expect_equal(res$y, c(3, 5) * 0.5)
})

test_that("runEmbedding (images) returns cells_id-ordered coordinate pairs", {
  # coords are out of cells_id order to exercise the ordering logic.
  # Seurat accessors are stubbed; meta.data comes from a real Seurat object.
  scdata <- mock_scdata()
  cells <- colnames(scdata)[1:3]
  # assign known, shuffled cells_id to the first 3 cells
  scdata@meta.data[cells, "cells_id"] <- c(2, 0, 1)

  coords <- data.frame(
    x = c(10, 20, 30),
    y = c(11, 21, 31),
    cell = cells,
    row.names = cells
  )

  stub(runEmbedding, "Seurat::Images", function(...) "fov1")
  stub(runEmbedding, "get_tissue_coords", function(img, sc) coords)

  req <- list(body = list(type = "images", use_saved = TRUE, config = list()))
  res <- runEmbedding(req, scdata)

  # length is max(cells_id) + 1
  expect_equal(length(res), max(scdata@meta.data$cells_id) + 1L)

  # cell with cells_id 0 is cells[2] -> (20,21)
  expect_equal(res[[1]], c(20, 21))
  # cells_id 1 is cells[3] -> (30,31)
  expect_equal(res[[2]], c(30, 31))
  # cells_id 2 is cells[1] -> (10,11)
  expect_equal(res[[3]], c(10, 11))
  # cells without a coordinate are NULL
  expect_null(res[[4]])
})

test_that("runEmbedding (polygons) dispatches to get_polygon_coords", {
  scdata <- mock_scdata()
  cells <- colnames(scdata)[1:2]
  scdata@meta.data[cells, "cells_id"] <- c(1, 0)

  coords <- data.frame(
    x = c(5, 7),
    y = c(6, 8),
    cell = cells,
    row.names = cells
  )

  mock_poly <- mock(coords)
  stub(runEmbedding, "Seurat::Images", function(...) "fov1")
  stub(runEmbedding, "get_polygon_coords", mock_poly)

  req <- list(body = list(type = "polygons", use_saved = TRUE, config = list()))
  res <- runEmbedding(req, scdata)

  expect_equal(length(mock_poly), 1)
  # cells_id 0 -> cells[2] -> (7,8); cells_id 1 -> cells[1] -> (5,6)
  expect_equal(res[[1]], c(7, 8))
  expect_equal(res[[2]], c(5, 6))
})

test_that("empty named list encodes consistently across JSON encoders", {
  # Test that empty named list (structure(list(), names = character(0)))
  # encodes identically in both jsonlite and yyjsonr for backward compatibility

  empty_named_list <- structure(list(), names = character(0))
  coordinate_pair <- c(1.5, 2.5)

  # Create test data with mixed content (coordinates and empty lists)
  test_data <- list(
    empty_named_list,  # Missing embedding
    coordinate_pair,   # Present embedding
    empty_named_list,  # Missing embedding
    c(3.1, 4.2)       # Present embedding
  )

  # Encode with jsonlite
  json_jsonlite <- as.character(jsonlite::toJSON(test_data, auto_unbox = TRUE))

  # Encode with yyjsonr
  json_yyjsonr <- yyjsonr::write_json_str(test_data, opts = list(
    dataframe = "columns",
    digits = 4,
    auto_unbox = TRUE
  ))

  # Both should produce identical JSON output
  expect_equal(json_jsonlite, json_yyjsonr)

  # Both should parse back identically
  parsed_jsonlite <- jsonlite::fromJSON(json_jsonlite)
  parsed_yyjsonr <- jsonlite::fromJSON(json_yyjsonr)

  expect_equal(parsed_jsonlite, parsed_yyjsonr)

  # Verify structure: empty entries should be empty lists, not NULL
  expect_equal(length(parsed_jsonlite[[1]]), 0)
  expect_equal(length(parsed_jsonlite[[3]]), 0)
  expect_true(is.list(parsed_jsonlite[[1]]))
  expect_true(is.list(parsed_jsonlite[[3]]))

  # Verify coordinates are preserved
  expect_equal(unname(parsed_jsonlite[[2]]), c(1.5, 2.5))
  expect_equal(unname(parsed_jsonlite[[4]]), c(3.1, 4.2))
})
