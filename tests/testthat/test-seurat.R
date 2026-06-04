# Seurat interoperability (skipped when Seurat is not installed).

test_that("as_sce converts a Seurat object to a SingleCellExperiment", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(1)
  m <- matrix(rpois(20 * 12, 5), nrow = 20,
              dimnames = list(paste0("g", 1:20), paste0("c", 1:12)))
  sce <- suppressWarnings({
    so <- SeuratObject::CreateSeuratObject(counts = m)
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    as_sce(so)
  })
  expect_true(is_sce(sce))
  expect_equal(ncol(sce), 12L)
  expect_true("counts" %in% SummarizedExperiment::assayNames(sce))
})

test_that("compute_qc runs on a Seurat-derived SCE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  set.seed(1)
  genes <- c("MT-CO1", paste0("GENE", 1:19))
  m <- matrix(rpois(20 * 12, 5), nrow = 20,
              dimnames = list(genes, paste0("c", 1:12)))
  sce <- suppressWarnings({
    so <- SeuratObject::CreateSeuratObject(counts = m)
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    as_sce(so)
  })

  qc <- compute_qc(sce, run_doublets = FALSE)
  cd <- as.data.frame(SummarizedExperiment::colData(qc))
  expect_true(all(c("nUMI", "nGene", "pct_mito") %in% names(cd)))
  expect_true(all(cd$pct_mito >= 0 & cd$pct_mito <= 100))
})
