# Render the reproducible QC report. The report's helper functions and data are
# injected into the knit environment (rather than referenced as cellboard::...),
# so the same code works whether cellboard is an installed package or has been
# sourced into the global environment by the Connect entry point (app.R).

#' Render a reproducible QC report
#'
#' Knits the bundled R Markdown template into a self-contained HTML report
#' describing the QC metrics, the chosen thresholds, and the retained-versus-
#' discarded outcome.
#'
#' @param sce An SCE. If it lacks QC metrics they are computed on the fly.
#' @param thresholds A `cb_thresholds` object or named list.
#' @param output_file Output file name (default `"cellboard_qc_report.html"`).
#' @param output_dir Directory to write into (default [tempdir()]).
#' @param title Report title.
#' @param sample_name Optional dataset/sample label shown in the report.
#' @param template Path to an `.Rmd` template (defaults to the bundled one).
#' @param quiet Passed to [rmarkdown::render()].
#'
#' @return The path to the rendered report.
#' @export
#'
#' @examples
#' \dontrun{
#'   sce <- compute_qc(example_sce())
#'   render_qc_report(sce, qc_thresholds(min_genes = 200, max_mito = 10))
#' }
render_qc_report <- function(sce,
                             thresholds = default_thresholds(),
                             output_file = "cellboard_qc_report.html",
                             output_dir = tempdir(),
                             title = "cellboard QC report",
                             sample_name = NULL,
                             template = NULL,
                             quiet = TRUE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Rendering the report requires the 'rmarkdown' package.", call. = FALSE)
  }
  sce <- ensure_counts(sce)
  th <- as_thresholds(thresholds)
  template <- template %||% cb_system_file("report", "qc_report.Rmd")
  if (!nzchar(template) || !file.exists(template)) {
    stop("Report template not found. Looked for inst/report/qc_report.Rmd.", call. = FALSE)
  }

  # Work in an isolated temp dir so we never write into a read-only package lib.
  work <- tempfile("cellboard_report_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  tmpl_local <- file.path(work, "qc_report.Rmd")
  file.copy(template, tmpl_local, overwrite = TRUE)

  # Inject data + helper functions into the knit environment.
  e <- new.env(parent = globalenv())
  e$sce             <- sce
  e$thresholds      <- th
  e$report_title    <- title
  e$sample_name     <- sample_name %||% "dataset"
  e$compute_qc      <- compute_qc
  e$apply_qc_filter <- apply_qc_filter
  e$summarise_qc    <- summarise_qc
  e$qc_metric_cols  <- qc_metric_cols
  e$qc_violin       <- qc_violin
  e$qc_scatter      <- qc_scatter
  e$qc_knee         <- qc_knee
  e$fmt_int         <- fmt_int

  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out <- rmarkdown::render(
    input       = tmpl_local,
    output_file = output_file,
    output_dir  = output_dir,
    envir       = e,
    quiet       = quiet
  )
  normalizePath(out, mustWork = FALSE)
}
