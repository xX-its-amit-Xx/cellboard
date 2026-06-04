# QC metric computation on a deterministic SCE.

test_that("compute_qc produces exact metric values", {
  sce <- make_test_sce()
  qc <- compute_qc(sce, run_doublets = FALSE)
  cd <- as.data.frame(SummarizedExperiment::colData(qc))

  expect_equal(cd$nUMI, c(10, 20, 20))
  expect_equal(cd$nGene, c(2L, 2L, 2L))
  expect_equal(cd$pct_mito, c(50, 0, 50))
})

test_that("mito detection uses the Symbol rowData column", {
  sce <- make_test_sce()
  qc <- compute_qc(sce, run_doublets = FALSE)
  meta <- S4Vectors::metadata(qc)$cellboard$qc
  expect_equal(meta$n_mito_genes, 1L)
  expect_equal(meta$mito_pattern, "^[Mm][Tt]-")
})

test_that("compute_qc warns and sets pct_mito = 0 when no mito genes match", {
  sce <- make_test_sce()
  expect_warning(
    qc <- compute_qc(sce, mito_pattern = "^NOMATCH-", run_doublets = FALSE),
    "No mitochondrial genes"
  )
  expect_true(all(as.data.frame(SummarizedExperiment::colData(qc))$pct_mito == 0))
})

test_that("compute_qc sets NA doublet metrics when run_doublets = FALSE", {
  qc <- compute_qc(make_test_sce(), run_doublets = FALSE)
  cd <- as.data.frame(SummarizedExperiment::colData(qc))
  expect_true(all(is.na(cd$doublet_score)))
  expect_true(all(is.na(cd$doublet_class)))
})

test_that("apply_qc_filter annotates and subsets the SCE", {
  qc <- compute_qc(make_test_sce(), run_doublets = FALSE)
  res <- apply_qc_filter(qc, qc_thresholds(max_mito = 10))
  expect_true("discard" %in% colnames(SummarizedExperiment::colData(res$sce)))
  # cells 1 and 3 have 50% mito -> discarded; cell 2 kept
  expect_equal(res$n_keep, 1L)
  expect_equal(ncol(res$filtered), 1L)
  expect_equal(res$keep, c(FALSE, TRUE, FALSE))
})

test_that("summarise_qc totals are internally consistent", {
  qc <- compute_qc(make_test_sce(), run_doublets = FALSE)
  res <- apply_qc_filter(qc, qc_thresholds(max_mito = 10))
  s <- summarise_qc(res)
  totals <- stats::setNames(s$totals$value, s$totals$metric)
  expect_equal(totals[["Retained"]] + totals[["Discarded"]], totals[["Input cells"]])
})

test_that("compute_qc works on a sparse counts matrix", {
  skip_if_not_installed("Matrix")
  sce <- make_test_sce()
  SummarizedExperiment::assay(sce, "counts") <-
    as(SummarizedExperiment::assay(sce, "counts"), "CsparseMatrix")
  qc <- compute_qc(sce, run_doublets = FALSE)
  expect_equal(as.data.frame(SummarizedExperiment::colData(qc))$nUMI,
               c(10, 20, 20))
})

test_that("scDblFinder runs end-to-end when requested", {
  skip_on_cran()
  skip_if_not_installed("scDblFinder")
  # needs a non-trivial number of cells; use the bundled example
  skip_if(!nzchar(system.file("extdata", "pbmc_small_sce.rds", package = "cellboard")))
  qc <- compute_qc(example_sce(), run_doublets = TRUE)
  cd <- as.data.frame(SummarizedExperiment::colData(qc))
  expect_true(is.numeric(cd$doublet_score))
  expect_true(any(!is.na(cd$doublet_score)))
})
