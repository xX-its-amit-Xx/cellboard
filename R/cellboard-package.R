#' cellboard: Interactive Single-Cell RNA-seq Quality Control
#'
#' cellboard is a Shiny dashboard for code-free quality control of single-cell
#' RNA-seq data. It loads 10x Genomics, [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#' or Seurat input, computes per-cell QC metrics (mitochondrial percentage,
#' genes and UMIs detected, and `scDblFinder` doublet scores), lets the user
#' set filtering thresholds interactively with live retained-versus-discarded
#' feedback, and exports a filtered object plus a reproducible report.
#'
#' The package is usable in three ways:
#' * launch the dashboard with [run_app()];
#' * call the QC engine programmatically ([compute_qc()], [apply_qc_filter()]);
#' * render a standalone report with [render_qc_report()].
#'
#' @keywords internal
"_PACKAGE"

# Standardised colData column names written by compute_qc() and consumed by the
# plotting helpers, the Shiny app and the report. Keeping them in one place means
# the whole pipeline agrees on what a metric is called. See qc_metric_cols().
.cb_cols <- c(
  nUMI          = "nUMI",
  nGene         = "nGene",
  pct_mito      = "pct_mito",
  doublet_score = "doublet_score",
  doublet_class = "doublet_class"
)

# Human-readable labels for the metrics, used in plots and the UI.
.cb_labels <- c(
  nUMI          = "UMIs per cell",
  nGene         = "Genes per cell",
  pct_mito      = "Mitochondrial %",
  doublet_score = "Doublet score",
  doublet_class = "Doublet class"
)

# These are column names referenced inside ggplot2::aes() (non-standard
# evaluation) in the plot helpers, not undefined globals. Declaring them keeps
# R CMD check quiet.
utils::globalVariables(c("group", "value", "status", "x", "y", "col",
                         "rank", "umi"))
