# A tiny, fully deterministic SingleCellExperiment for exact-value assertions.
#
# counts (genes x cells):
#            c1  c2  c3
#   MT-CO1    5   0  10
#   GENE1     5  10   0
#   GENE2     0  10   0
#   GENE3     0   0  10
#   GENE4/5   0   0   0
#
# => nUMI  = c(10, 20, 20)
#    nGene = c( 2,  2,  2)
#    pct_mito = c(50, 0, 50)

# Convenience predicate for tests (mirrors the internal cellboard helper).
is_sce <- function(x) methods::is(x, "SingleCellExperiment")

make_test_sce <- function() {
  testthat::skip_if_not_installed("SingleCellExperiment")
  m <- matrix(
    c(5, 0, 10,
      5, 10, 0,
      0, 10, 0,
      0, 0, 10,
      0, 0, 0,
      0, 0, 0),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:6), paste0("c", 1:3))
  )
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = m),
    rowData = S4Vectors::DataFrame(
      Symbol = c("MT-CO1", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
    )
  )
}
