stub_updateCellSetsThroughApi <- function(cell_sets_object,
                                          api_url,
                                          experiment_id,
                                          cell_set_key,
                                          auth_JWT,
                                          append = TRUE) {
  # empty function to simplify mocking. we test patching independently.
}

stubbed_cellCycleScoring <- function(req, data) {
  mockery::stub(cellCycleScoring,
                "updateCellSetsThroughApi",
                stub_updateCellSetsThroughApi)
  cellCycleScoring(req, data)
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
    new_names[41:121] <- Seurat::cc.genes$g2m.genes[1:80]
  }

  gene_annotations <- data.frame(
    input = enids,
    name = new_names,
    original_name = row.names(counts),
    row.names = enids
  )

  pbmc_small <- SeuratObject::CreateSeuratObject(counts = counts)

  if (phase == "S") {
    pbmc_small@assays$RNA@counts[1:40, ] <- 3
    pbmc_small@assays$RNA@counts[41:121, ] <- 0
  } else{
    pbmc_small@assays$RNA@counts[41:121, ] <- 3
  }


  pbmc_small$cells_id <- 0:(ncol(pbmc_small) - 1)
  pbmc_small@misc$gene_annotations <- gene_annotations
  pbmc_small@misc$color_pool <- mock_color_pool(20)

  pbmc_small <-
    Seurat::NormalizeData(pbmc_small,
                          normalization.method = "LogNormalize",
                          verbose = FALSE)
  pbmc_small <-
    Seurat::FindVariableFeatures(pbmc_small, verbose = FALSE)
  pbmc_small <- Seurat::ScaleData(pbmc_small, verbose = FALSE)
  return(pbmc_small)
}

mock_req <- function(data) {
  req <- list(body = list())

  return(req)
}

test_that("Cell Cycle Scoring returns expected and formatted result", {
  scdata <- mock_scdata()
  req <- mock_req()

  result <- stubbed_cellCycleScoring(scdata, req)

  expect_true(result$key == "Phase")
  expect_true(result$name == "Phase")
  expect_true(result$type == "cellSets")

  expect_equal(length(result$children), 1)
  expect_true(result$children[[1]]$key == "Phase-S")
  expect_true(result$children[[1]]$name == "S")
})

# Start writing tests
test_that("run_cell_cycle_scoring returns a data frame with correct columns",
          {
            scdata <- mock_scdata()

            result <- run_cell_cycle_scoring(scdata, req)

            expect_is(result, "data.frame")

            expect_true("cluster" %in% colnames(result))
            expect_true("cell_ids" %in% colnames(result))
          })

test_that("run_cell_cycle_scoring properly classifies S cells", {
  scdata <- mock_scdata("S")

  result <- run_cell_cycle_scoring(scdata, req)

  expect_is(result, "data.frame")

  expect_true(all(result$cluster == "S"))
})

test_that("run_cell_cycle_scoring properly classifies G2M cells", {
  scdata <- mock_scdata("G2M")

  result <- run_cell_cycle_scoring(scdata, req)

  expect_is(result, "data.frame")

  expect_true(all(result$cluster == "G2M"))
})

test_that("run_cell_cycle_scoring returns 'Undetermined' cellset when no genes are detected.",
          {
            mock_scdata <- mock_scdata(phase = "S", add_cycle_genes = FALSE)

            result <- run_cell_cycle_scoring(mock_scdata, req)

            expect_is(result, "data.frame")

            expect_true(all(result$cluster == "Undetermined"))
          })

test_that("format_cluster_cellsets returns a list with the correct structure",
          {
            # Create a mock data frame for testing
            mock_cell_sets <- data.frame(cluster = c("A", "A", "B", "C", "C"),
                                         cell_ids = 1:5)

            # Call the function
            result <-
              format_cluster_cellsets(mock_cell_sets, "method", c("red", "blue", "green"))

            # Check if the result is a list
            expect_is(result, "list")

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
            expect_is(result$children, "list")

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

test_that("format_cluster_cellsets handles various cases", {
  mock_cell_sets <- data.frame(cluster = c("A", "A", "B", "C", "C"),
                               cell_ids = 1:5)

  result <-
    format_cluster_cellsets(mock_cell_sets, "method", c("red", "blue", "green"))

  # Key Values
  expect_equal(result$key, "method")
  expect_equal(result$children[[1]]$key, "method-A")

  # Color and cell id Assignment
  expect_equal(result$children[[1]]$color, "red")
  expect_equal(result$children[[2]]$color, "blue")
  expect_equal(result$children[[1]]$cellIds, c(1, 2))
  expect_equal(result$children[[2]]$cellIds, list(3))

  # Numeric Cluster Labels
  mock_cell_sets$cluster <- c(1, 1, 2, 3, 3)
  result_numeric <-
    format_cluster_cellsets(mock_cell_sets, "method", c("red", "blue", "green"))
  expect_equal(result_numeric$children[[1]]$name, "Cluster 1")
})
