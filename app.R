# app.R — entry point for shiny::runApp() and Posit Connect.
#
# The deployment bundle ships the R/ sources, so we load those directly: a fresh
# clone (or a Connect bundle) runs without first installing the cellboard
# package. Dependencies (shiny, bslib, Bioconductor, ...) are restored from
# manifest.json / renv.lock. If the sources are absent but cellboard is
# installed, fall back to the installed package.

local({
  r_dir <- "R"
  r_files <- if (dir.exists(r_dir)) {
    list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  } else character(0)

  if (length(r_files)) {
    invisible(lapply(sort(r_files), source))
  } else if (requireNamespace("cellboard", quietly = TRUE)) {
    library(cellboard)
  } else {
    stop("cellboard sources (R/) not found and the package is not installed.")
  }
})

# The trailing shiny.appobj is what runApp()/Connect serve.
run_app()
