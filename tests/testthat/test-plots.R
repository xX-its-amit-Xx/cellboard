# Plot helpers return ggplot objects and tolerate missing metrics.

test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("ggplot2")
  qc <- compute_qc(example_sce(), run_doublets = FALSE)
  res <- apply_qc_filter(qc, default_thresholds(qc))
  th <- default_thresholds(qc)

  expect_s3_class(qc_violin(res$sce, "pct_mito", th), "ggplot")
  expect_s3_class(qc_violin(res$sce, "nUMI", th), "ggplot")
  expect_s3_class(qc_scatter(res$sce, "nUMI", "nGene", "status", th), "ggplot")
  expect_s3_class(qc_knee(res$sce, th), "ggplot")
})

test_that("qc_violin handles an all-NA metric gracefully", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("ggplot2")
  qc <- compute_qc(example_sce(), run_doublets = FALSE) # doublet_score all NA
  expect_s3_class(qc_violin(qc, "doublet_score"), "ggplot")
})
