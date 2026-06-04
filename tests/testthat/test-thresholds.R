# Pure threshold / filtering logic. These run without Bioconductor.

test_that("qc_thresholds validates inputs", {
  th <- qc_thresholds(min_genes = 200, max_mito = 10)
  expect_s3_class(th, "cb_thresholds")
  expect_equal(th$min_genes, 200)
  expect_null(th$max_genes)

  expect_error(qc_thresholds(min_genes = "x"), "single number")
  expect_error(qc_thresholds(min_genes = c(1, 2)), "single number")
  expect_error(qc_thresholds(min_genes = 500, max_genes = 100), "min_genes")
  expect_error(qc_thresholds(min_umi = 1000, max_umi = 10), "min_umi")
})

test_that("default_thresholds returns sensible defaults", {
  th <- default_thresholds()
  expect_s3_class(th, "cb_thresholds")
  expect_true(th$min_genes > 0)
  expect_true(th$max_mito > 0 && th$max_mito <= 100)
  expect_true(isTRUE(th$remove_doublet_class))
})

test_that("filter_cells_df keeps everything with empty thresholds", {
  df <- data.frame(nUMI = c(100, 5000), nGene = c(50, 1500),
                   pct_mito = c(30, 5))
  res <- filter_cells_df(df, qc_thresholds())
  expect_true(all(res$keep))
  expect_equal(ncol(res$reasons), 0L)
})

test_that("filter_cells_df applies each threshold", {
  df <- data.frame(
    nUMI          = c(100, 5000, 8000, 3000),
    nGene         = c(40,  1500, 2200, 800),
    pct_mito      = c(60,  4,    5,    50),
    doublet_class = c("singlet", "singlet", "doublet", "singlet")
  )
  th <- qc_thresholds(min_genes = 200, min_umi = 500, max_mito = 10,
                      remove_doublet_class = TRUE)
  res <- filter_cells_df(df, th)
  # row1 fails genes+umi+mito; row2 passes; row3 fails doublet; row4 fails mito
  expect_equal(res$keep, c(FALSE, TRUE, FALSE, FALSE))
  expect_true("high_mito" %in% names(res$reasons))
  expect_true("doublet_class" %in% names(res$reasons))
})

test_that("NA metrics never cause a discard", {
  df <- data.frame(nUMI = c(NA, 5000), nGene = c(NA, 1500),
                   pct_mito = c(NA, 5))
  res <- filter_cells_df(df, qc_thresholds(min_genes = 200, max_mito = 10))
  expect_true(all(res$keep))
})

test_that("stricter thresholds are monotone (retain no more cells)", {
  df <- data.frame(nUMI = c(100, 1000, 5000), nGene = c(50, 500, 2000),
                   pct_mito = c(40, 8, 3))
  loose  <- filter_cells_df(df, qc_thresholds(min_genes = 100))
  strict <- filter_cells_df(df, qc_thresholds(min_genes = 600, max_mito = 5))
  expect_lte(sum(strict$keep), sum(loose$keep))
})

test_that("filter_cells_df accepts a plain list of thresholds", {
  df <- data.frame(nUMI = 100, nGene = 50, pct_mito = 60)
  res <- filter_cells_df(df, list(max_mito = 10))
  expect_false(res$keep)
})
