# Input loading and type detection.

test_that("example_sce returns a valid SCE with counts", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- example_sce()
  expect_true(is_sce(sce))
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
  expect_gt(ncol(sce), 0)
  expect_gt(nrow(sce), 0)
})

test_that("detect_input_type recognises the bundled inputs", {
  skip_if_not_installed("SingleCellExperiment")
  dir10x <- example_sce(as = "10x_path")
  expect_equal(detect_input_type(dir10x), "10x_dir")

  rds <- system.file("extdata", "pbmc_small_sce.rds", package = "cellboard")
  skip_if(!nzchar(rds))
  expect_equal(detect_input_type(rds), "rds")
})

test_that("detect_input_type errors on unsupported input", {
  expect_error(detect_input_type(""), "non-empty")
  tmp <- tempfile(fileext = ".txt")
  writeLines("hello", tmp)
  expect_error(detect_input_type(tmp), "Unsupported")
})

test_that("read_10x round-trips the bundled 10x directory", {
  skip_if_not_installed("DropletUtils")
  sce_rds <- example_sce()
  sce_10x <- read_10x(example_sce(as = "10x_path"))
  expect_equal(ncol(sce_10x), ncol(sce_rds))
  expect_equal(nrow(sce_10x), nrow(sce_rds))
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce_10x))
})

test_that("read_input dispatches on type and validates", {
  skip_if_not_installed("DropletUtils")
  sce <- read_input(example_sce(as = "10x_path"))      # auto -> 10x
  expect_true(is_sce(sce))
})

test_that("read_input descends into a nested 10x matrix directory", {
  skip_if_not_installed("DropletUtils")
  src <- example_sce(as = "10x_path")                  # flat 10x dir
  parent <- tempfile("nested10x_")
  leaf <- file.path(parent, "hg19")
  dir.create(leaf, recursive = TRUE)
  file.copy(list.files(src, full.names = TRUE), leaf)
  # files live in parent/hg19/, but the user points at the parent
  expect_equal(detect_input_type(parent), "10x_dir")
  sce <- read_input(parent)
  expect_true(is_sce(sce))
  expect_gt(ncol(sce), 0)
})

test_that("as_sce coerces a bare matrix", {
  skip_if_not_installed("SingleCellExperiment")
  m <- matrix(1:12, nrow = 3)
  sce <- as_sce(m)
  expect_true(is_sce(sce))
  expect_equal(dim(sce), c(3L, 4L))
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
})

test_that("ensure_counts renames a sole assay to 'counts'", {
  skip_if_not_installed("SingleCellExperiment")
  m <- matrix(1:6, nrow = 2)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(foo = m))
  out <- as_sce(sce)
  expect_true("counts" %in% SummarizedExperiment::assayNames(out))
})

test_that("read_seurat errors clearly on a non-Seurat file", {
  skip_if_not_installed("SingleCellExperiment")
  tmp <- tempfile(fileext = ".rds")
  saveRDS(example_sce(), tmp)
  expect_error(read_seurat(tmp), "does not contain a Seurat")
})
