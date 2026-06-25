# IMPORTANT: these tests are duplicated in the pipeline. If you update, change both.

mock_req <- function(type = "louvain") {
  req <- list(body =
                list(
                  config = list(
                    resolution = 2,
                    apiUrl = "mock_api_url",
                    authJwt = "mock_auth"

                  ),
                  type = type
                ))
}


mock_cellset_object <- function(n_cells, n_clusters) {
  # cellset objects are data.frames. They are formatted to lists by
  # format_cell_sets_object

  if (n_clusters == 0) {
    return(data.frame(cluster = integer(), cell_ids = integer()))
  }

  # cell_ids with no replacement to avoid repeated cell_ids
  data.frame(cluster = sample(1:n_clusters, size = n_cells, replace = T),
             cell_ids = sample(1:(2*n_cells), size = n_cells))
}


stub_updateCellSetsThroughApi <- function(cell_sets_object,
                                          api_url,
                                          experiment_id,
                                          cell_set_key,
                                          auth_JWT,
                                          append = TRUE) {

  # empty function to simplify mocking. we test patching independently.
}

stubbed_runClusters <- function(req, data) {
  mockery::stub(runClusters,
                "updateCellSetsThroughApi",
                stub_updateCellSetsThroughApi)
  runClusters(req, data)
}

test_that("runClusters works with bpcells", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata(use_bpcells = TRUE)
  expected_keys <- c("cluster", "cell_ids")

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, data)
    expect_equal(names(res), expected_keys)
  }
})

test_that("runClusters returns correct keys", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_keys <- c("cluster", "cell_ids")

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, data)
    expect_equal(names(res), expected_keys)
  }
})

test_that("runClusters returns one value per cell", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_n_cells <- ncol(data)

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, data)
    n_cells <- nrow(res)

    expect_equal(n_cells, expected_n_cells)
  }
})

test_that("runClusters orders barcodes correctly", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()
  expected_barcodes <- colnames(data)

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, data)
    barcodes <- rownames(res)
    expect_equal(barcodes, expected_barcodes)
  }
})

test_that("runClusters returns at least one cluster", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, data)
    n_clusters <- length(unique(res$cluster))
    expect_gte(n_clusters, 1)
  }
})

test_that("runClusters uses active.reduction in misc slot", {
  algos <- c("louvain", "leiden")
  data <- mock_scdata()

  # get error if no PCA/SNN graph
  blah_reduction <- data@reductions$pca
  data@reductions$pca <- NULL
  data@graphs$RNA_snn <- NULL

  for (algo in algos) {
    req <- mock_req(type = algo)
    expect_error(stubbed_runClusters(req, data), "'pca' not found")
  }

  # will use active.reduction to get SNN graph
  data@reductions$blah_reduction <- blah_reduction
  data@misc$active.reduction <- "blah_reduction"

  for (algo in algos) {
    req <- mock_req(type = algo)
    expect_error(stubbed_runClusters(req, data), NA)
  }
})

with_fake_http(test_that("updateCellSetsThroughApi sends patch request", {
  expect_PATCH(
    updateCellSetsThroughApi(list(), "api_url", "experiment_id", "cell_set_key", "auth")
  )
}))


test_that("format_cell_sets_object returns correct items", {
  n_clusters <- 5
  cell_sets <- mock_cellset_object(1000, n_clusters)
  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  types <- c(
    "key" = "character",
    "name" = "character",
    "rootNode" = "logical",
    "type" = "character",
    "children" = "list"
  )

  for (algo in algos) {
    res <- format_cell_sets_object(cell_sets, algo, color_pool)

    # expect all keys present
    expect_setequal(names(res), names(types))

    for (item in names(types)) {
      expect_type(res[[item]], types[[item]])
    }
  }
})


test_that("format_cell_sets_object correctly formats a cellset object", {

  n_clusters <- 5
  cell_sets <- mock_cellset_object(1000, n_clusters)
  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  expected_items <- c("key", "name", "rootNode", "type", "color", "cellIds")

  for (algo in algos) {
    res <- format_cell_sets_object(cell_sets, algo, color_pool)


    # each children has the expected items
    for (cluster in res$children) {
      expect_setequal(names(cluster), expected_items)
    }

    # each cluster has correct content
    for (cluster in res$children) {
      expect_true(startsWith(cluster$key, algo))
      expect_true(startsWith(cluster$name, "Cluster"))
      expect_type(cluster$cellIds, "integer")
      expect_gte(length(cluster$cellIds), 1)
    }
  }
})


test_that("format_cell_sets_object result has correct number of clusters",{
  n_clusters <- c(0, 1, 4, 6)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    for (n in n_clusters) {
    cell_sets <- mock_cellset_object(100, n)
    color_pool <- mock_color_pool(n)

    res <- format_cell_sets_object(cell_sets, algo, color_pool)

    # number of clusters is the number of children elements
    expect_equal(length(res$children), n)

  }}
})

test_that("format_cell_sets_object returns empty children on empty cellset", {
  n_clusters <- 5
  cell_sets <- mock_cellset_object(0, n_clusters)

  color_pool <- mock_color_pool(n_clusters)
  algos <- c("louvain", "leiden")

  for (algo in algos) {
    res <- format_cell_sets_object(cell_sets, algo, color_pool)
    expect_equal(length(res$children), 0)
  }
})


test_that("format_cell_sets_object consumes palette in sorted cluster order", {
  # the palette (e.g. from Spaco) is ordered to match sort(unique(cluster))
  cell_sets <- data.frame(
    cluster = c(2, 2, 1, 3),
    cell_ids = c(10, 11, 20, 30)
  )
  palette <- c("#aaaaaa", "#bbbbbb", "#cccccc")

  res <- format_cell_sets_object(cell_sets, "louvain", palette)
  assigned <- vapply(res$children, function(x) x$color, character(1))

  # clusters iterate as 1, 2, 3 and take palette[1], palette[2], palette[3]
  expect_equal(assigned, palette)
})


test_that("get_spaco_color_map returns NULL when there are no slices", {
  # no images -> returns before any distance/embedding work
  data <- mock_scdata()
  cell_sets <- mock_cellset_object(100, 4)

  mockery::stub(
    get_spaco_color_map, "spatial_distance_r",
    function(...) stop("should not compute distances without slices")
  )
  expect_null(get_spaco_color_map(data, cell_sets))
})


mock_spatial_grid <- function(side, n_clusters) {
  coords <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  labels <- as.character(rep(seq_len(n_clusters), length.out = nrow(coords)))
  list(coords = coords, labels = labels)
}


test_that("spatial_distance_r returns a symmetric zero-diagonal matrix", {
  grid <- mock_spatial_grid(10, 4)
  m <- spatial_distance_r(grid$coords, grid$labels, radius = 5)

  expect_equal(dim(m), c(4, 4))
  expect_equal(rownames(m), as.character(1:4))
  expect_true(isSymmetric(unname(m)))
  expect_equal(unname(diag(m)), rep(0, 4))
})


test_that("spatial_distance_r matches Spaco's interlacement on a grid", {
  # spatial_distance_r is a port of spaco.distance.spatial_distance; check the
  # score for a known layout (cluster A and B interlaced in a checkerboard)
  coords <- as.matrix(expand.grid(x = 1:6, y = 1:6))
  labels <- ifelse((coords[, 1] + coords[, 2]) %% 2 == 0, "A", "B")
  m <- spatial_distance_r(coords, labels, radius = 1.5, n_cells = 1L)

  # checkerboard: every neighbour within radius 1.5 is the other cluster, so
  # A-B interlacement is non-zero and the diagonal is zero
  expect_gt(m["A", "B"], 0)
  expect_equal(m["A", "B"], m["B", "A"])
  expect_equal(unname(diag(m)), c(0, 0))
})


test_that("embed_graph_r returns one hex color per cluster", {
  grid <- mock_spatial_grid(12, 6)
  m <- spatial_distance_r(grid$coords, grid$labels, radius = 4)
  cols <- embed_graph_r(m)

  expect_named(cols, rownames(m))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
  expect_equal(length(unique(cols)), length(cols))
})


test_that("embed_graph_r errors for too few clusters to embed", {
  # with 2 clusters n_neighbors collapses to 1 and uwot's spectral init errors;
  # get_spaco_color_map's tryCatch turns this into a fall back to the color pool
  grid <- mock_spatial_grid(12, 2)
  m <- spatial_distance_r(grid$coords, grid$labels, radius = 4)
  expect_error(embed_graph_r(m))
})


test_that("merge_cluster_distances sums slices aligned to the cluster set", {
  m1 <- matrix(c(0, 1, 1, 0), 2, dimnames = list(c("1", "2"), c("1", "2")))
  m2 <- matrix(c(0, 2, 2, 0), 2, dimnames = list(c("2", "3"), c("2", "3")))

  merged <- merge_cluster_distances(list(m1, m2), c("1", "2", "3"))

  expect_equal(dim(merged), c(3, 3))
  expect_equal(merged["1", "2"], 1)
  expect_equal(merged["2", "3"], 2)
  expect_equal(merged["1", "3"], 0)
})


# Build the (cell_sets, coords) pair get_spaco_color_map consumes from a grid:
# cell_sets is the data.frame indexed by cell position, coords is what a stubbed
# GetTissueCoordinates returns (a cell/x/y frame in the same order).
mock_spaco_inputs <- function(side, n_clusters) {
  grid <- mock_spatial_grid(side, n_clusters)
  n <- nrow(grid$coords)
  list(
    cell_sets = data.frame(cluster = as.integer(grid$labels), cell_ids = seq_len(n)),
    coords = data.frame(cell = seq_len(n), x = grid$coords[, "x"], y = grid$coords[, "y"])
  )
}


test_that("get_spaco_color_map returns one hex colour per cluster on the happy path", {
  fixtures <- mock_spaco_inputs(12, 6)
  data <- mock_scdata()

  # stub the slice accessors so no real FOV/VisiumV2 object is needed; the
  # distance/embedding maths (spatial_distance_r -> embed_graph_r) runs for real
  mockery::stub(get_spaco_color_map, "Seurat::Images", function(...) "fov")
  mockery::stub(get_spaco_color_map, "get_img_scale", function(...) NULL)
  mockery::stub(
    get_spaco_color_map, "SeuratObject::GetTissueCoordinates",
    function(...) fixtures$coords
  )

  colors <- get_spaco_color_map(data, fixtures$cell_sets)

  expect_length(colors, 6)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})


test_that("get_spaco_color_map falls back to NULL when embedding fails", {
  # 2 clusters collapse uwot's spectral init -> embed_graph_r errors -> the
  # tryCatch returns NULL so runClusters keeps the default colour pool
  fixtures <- mock_spaco_inputs(12, 2)
  data <- mock_scdata()

  mockery::stub(get_spaco_color_map, "Seurat::Images", function(...) "fov")
  mockery::stub(get_spaco_color_map, "get_img_scale", function(...) NULL)
  mockery::stub(
    get_spaco_color_map, "SeuratObject::GetTissueCoordinates",
    function(...) fixtures$coords
  )

  expect_null(get_spaco_color_map(data, fixtures$cell_sets))
})


test_that("runClusters does not crash with less than 10 dimensions available", {
  algos <- c("louvain", "leiden")
  scdata <- mock_scdata()
  expected_keys <- c("cluster", "cell_ids")

  # remove all pre-existing reductions and calculate low-PC PCA
  scdata <- Seurat::DietSeurat(scdata)
  scdata <- suppressWarnings(Seurat::RunPCA(scdata, assay = "RNA", npcs = 2, verbose = F))

  for (algo in algos) {
    req <- mock_req(type = algo)
    res <- stubbed_runClusters(req, scdata)
    expect_equal(names(res), expected_keys)
  }
})


test_that("getClusters uses the default value of 10 if there are enough PCs available",{
  algos <- c("louvain", "leiden")
  scdata <- mock_scdata()
  resolution <- 0.8

  # remove all pre-existing reductions and calculate low-PC PCA
  scdata <- Seurat::DietSeurat(scdata)
  scdata@commands <- list()
  scdata <- suppressWarnings(Seurat::RunPCA(scdata, assay = "RNA", npcs = 20, verbose = F))

  for (algo in algos) {
    clustered_scdata <- getClusters(algo, resolution, scdata)
    if (algo == "louvain") expect_equal(clustered_scdata@commands$FindNeighbors.RNA.pca$dims, 1:10)
    # difficult to test in leiden, so test internal state as proxy
    if (algo == "leiden") expect_true("seurat_clusters" %in% names(clustered_scdata@meta.data))
  }
})
