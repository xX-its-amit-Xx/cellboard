#!/usr/bin/env Rscript
# Dev helper: install the cellboard runtime + test stack.
# Logs progress to dev/install_log.txt and writes dev/install_status.json at the end.
# Safe to re-run: already-installed packages are skipped.

options(warn = 1)
`%||%` <- function(a, b) if (is.null(a)) b else a

repo_cran <- "https://cloud.r-project.org"
options(repos = c(CRAN = repo_cran), timeout = 1200, Ncpus = max(1, parallel::detectCores() - 2))

log_path <- "dev/install_log.txt"
status_path <- "dev/install_status.json"

say <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = log_path, append = TRUE, sep = "")
}

installed_now <- function() rownames(installed.packages())

ensure_cran <- function(pkgs) {
  have <- installed_now()
  need <- setdiff(pkgs, have)
  if (!length(need)) { say("CRAN ok (already present): ", paste(pkgs, collapse = ", ")); return(invisible()) }
  say("Installing CRAN: ", paste(need, collapse = ", "))
  tryCatch(
    install.packages(need, quiet = TRUE),
    error = function(e) say("  ERROR installing CRAN [", paste(need, collapse = ","), "]: ", conditionMessage(e))
  )
}

ensure_bioc <- function(pkgs) {
  have <- installed_now()
  need <- setdiff(pkgs, have)
  if (!length(need)) { say("Bioc ok (already present): ", paste(pkgs, collapse = ", ")); return(invisible()) }
  say("Installing Bioc: ", paste(need, collapse = ", "))
  tryCatch(
    BiocManager::install(need, update = FALSE, ask = FALSE, quiet = TRUE),
    error = function(e) say("  ERROR installing Bioc [", paste(need, collapse = ","), "]: ", conditionMessage(e))
  )
}

unlink(log_path)
say("=== cellboard dependency install starting ===")
say("R: ", R.version.string, " | Ncpus=", getOption("Ncpus"))

# 1. BiocManager bootstrap
ensure_cran("BiocManager")
suppressMessages(suppressWarnings(require(BiocManager)))
if (requireNamespace("BiocManager", quietly = TRUE)) {
  say("Bioconductor version: ", as.character(BiocManager::version()))
}

# 2. Lightweight CRAN UI + utility deps (fast, makes the app shell runnable)
ensure_cran(c("shiny", "DT", "plotly", "shinycssloaders", "shinyWidgets", "bsicons", "jsonlite"))

# 3. Core Bioconductor single-cell stack (the QC engine)
ensure_bioc(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
              "scater", "scuttle", "scran", "scDblFinder", "DropletUtils"))

# 4. Optional heavy interop (best effort; cellboard treats Seurat as Suggests)
ensure_cran(c("SeuratObject", "Seurat"))

# Final status report
core <- c("shiny","bslib","DT","plotly","SingleCellExperiment","SummarizedExperiment",
          "scater","scuttle","scran","scDblFinder","DropletUtils","Matrix","testthat","rmarkdown")
optional <- c("SeuratObject","Seurat","shinycssloaders","shinyWidgets","bsicons")
have <- installed_now()
status <- list(
  finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  core_installed = as.list(setNames(core %in% have, core)),
  optional_installed = as.list(setNames(optional %in% have, optional)),
  core_all_present = all(core %in% have)
)
writeLines(jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE), status_path)
say("Core all present: ", all(core %in% have))
say("Missing core: ", paste(setdiff(core, have), collapse = ", "))
say("=== DONE ===")
