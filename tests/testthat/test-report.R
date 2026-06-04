# Reproducible report rendering (slow; skipped on CRAN).

test_that("render_qc_report produces a non-trivial HTML file", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("SingleCellExperiment")
  skip_if(!rmarkdown::pandoc_available(), "pandoc not available")
  skip_if(!nzchar(system.file("report", "qc_report.Rmd", package = "cellboard")))

  qc <- compute_qc(example_sce(), run_doublets = FALSE)
  out <- render_qc_report(
    qc, qc_thresholds(min_genes = 200, max_mito = 10),
    output_file = "test_report.html", output_dir = tempdir(), quiet = TRUE
  )
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 1000)
})
