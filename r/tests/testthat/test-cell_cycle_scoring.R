stub_updateCellSetsThroughApi <- function(cell_sets_object,
                                          api_url,
                                          experiment_id,
                                          cell_set_key,
                                          auth_JWT,
                                          append = TRUE) {
  # empty function to simplify mocking. we test patching independently.
}

stubbed_cellCycleScoring <- function(req, scdata) {
  mockery::stub(
    cellCycleScoring,
    "updateCellSetsThroughApi",
    stub_updateCellSetsThroughApi
  )
  suppressWarnings(cellCycleScoring(req, scdata))
}

mock_color_pool <- function(n) {
  paste0("color_", 1:n)
}

mock_scdata <- function(phase = "S",
                        add_cycle_genes = TRUE) {
  enids <- paste0("ENSG", 1:3000)
  cell_names <- paste0("cell_", 1:500)

  set.seed(ULTIMATE_SEED)
  counts <-
    matrix(
      rnbinom(
        length(enids) * length(cell_names),
        mu = 0.3,
        size = 5
      ),
      ncol = length(cell_names),
      byrow = TRUE
    )
  rownames(counts) <- enids
  colnames(counts) <- cell_names

  new_names <- enids
  if (add_cycle_genes) {
    new_names[1:40] <- Seurat::cc.genes$s.genes[1:40]
    new_names[41:90] <- Seurat::cc.genes$g2m.genes[1:50]
  }

  gene_annotations <- data.frame(
    input = enids,
    name = new_names,
    original_name = row.names(counts),
    row.names = enids
  )

  scdata <- SeuratObject::CreateSeuratObject(counts = counts)

  # We are forcing expression numbers to return expected cell cycle phase.
  if (phase == "S") {
    scdata@assays$RNA@counts[1:40, ] <- 3
    scdata@assays$RNA@counts[41:90, ] <- 0
  } else {
    scdata@assays$RNA@counts[41:90, ] <- 3
  }


  scdata$cells_id <- 0:(ncol(scdata) - 1)
  scdata@misc$gene_annotations <- gene_annotations
  scdata@misc$color_pool <- mock_color_pool(20)

  scdata <-
    Seurat::NormalizeData(scdata,
      normalization.method = "LogNormalize",
      verbose = FALSE
    )
  scdata <-
    Seurat::FindVariableFeatures(scdata, verbose = FALSE)
  scdata <- Seurat::ScaleData(scdata, verbose = FALSE)
  return(scdata)
}

mock_req <- function(data) {
  req <- list(body = list())

  return(req)
}

test_that("Cell Cycle Scoring returns expected and formatted result", {
  scdata <- mock_scdata()
  req <- mock_req()

  result <- stubbed_cellCycleScoring(req, scdata)

  expect_true(result$key == "Phase")
  expect_true(result$name == "Phase")
  expect_true(result$type == "cellSets")

  expect_equal(length(result$children), 1)
  expect_true(result$children[[1]]$key == "Phase-S")
  expect_true(result$children[[1]]$name == "S")
})

test_that("run_cell_cycle_scoring returns a data frame with correct columns", {
  scdata <- mock_scdata()

  result <- suppressWarnings(run_cell_cycle_scoring(scdata))

  expect_s3_class(result, "data.frame")

  expect_true("cluster" %in% colnames(result))
  expect_true("cell_ids" %in% colnames(result))
})

test_that("run_cell_cycle_scoring properly classifies S cells", {
  scdata <- mock_scdata("S")

  result <- suppressWarnings(run_cell_cycle_scoring(scdata))

  expect_s3_class(result, "data.frame")

  expect_true(all(result$cluster == "S"))
})

test_that("run_cell_cycle_scoring properly classifies G2M cells", {
  scdata <- mock_scdata("G2M")

  result <- suppressWarnings(run_cell_cycle_scoring(scdata))

  expect_s3_class(result, "data.frame")

  expect_true(all(result$cluster == "G2M"))
})

test_that("run_cell_cycle_scoring returns 'Undetermined' cellset when no genes are detected.", {
  scdata <- mock_scdata(phase = "S", add_cycle_genes = FALSE)

  result <- suppressWarnings(run_cell_cycle_scoring(scdata))

  expect_s3_class(result, "data.frame")

  expect_true(all(result$cluster == "Undetermined"))
})

test_that("format_phase_cellsets returns a list with the correct structure", {
  # Create a mock data frame for testing
  mock_cell_sets <- data.frame(
    cluster = c("A", "A", "B", "C", "C"),
    cell_ids = 1:5
  )

  # Call the function
  result <-
    format_cluster_cellsets(mock_cell_sets, c("red", "blue", "green"))

  # Check if the result is a list
  expect_type(result, "list")

  # Check if the list has the correct elements
  expect_true("key" %in% names(result))
  expect_true("name" %in% names(result))
  expect_true("rootNode" %in% names(result))
  expect_true("type" %in% names(result))
  expect_true("children" %in% names(result))

  # Check the rootNode value
  expect_true(result$rootNode)

  # Check the type value
  expect_equal(result$type, "cellSets")

  # Check if the children element is a list
  expect_type(result$children, "list")

  # Check the structure of each child
  for (child in result$children) {
    expect_true("key" %in% names(child))
    expect_true("name" %in% names(child))
    expect_true("rootNode" %in% names(child))
    expect_true("type" %in% names(child))
    expect_true("color" %in% names(child))
    expect_true("cellIds" %in% names(child))
  }
})

test_that("format_phase_cellsets returns expected formatting", {
  mock_cell_sets <- data.frame(
    cluster = c("A", "A", "B", "C", "C"),
    cell_ids = 1:5
  )

  result <-
    format_phase_cellsets(mock_cell_sets, c("red", "blue", "green"))

  # Key Values
  expect_equal(result$key, "Phase")
  expect_equal(result$children[[1]]$key, "Phase-A")

  # Color and cell id Assignment
  expect_equal(result$children[[1]]$color, "red")
  expect_equal(result$children[[2]]$color, "blue")
  expect_equal(result$children[[1]]$cellIds, c(1, 2))
  expect_equal(result$children[[2]]$cellIds, list(3))
})
