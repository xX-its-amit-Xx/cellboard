#!/usr/bin/env Rscript
# Generate the Posit Connect deployment artifacts:
#   manifest.json  - via rsconnect::writeManifest() (scans app.R + R/)
#   renv.lock      - runtime dependencies pinned to the installed versions
#
# Re-run this whenever the app's runtime dependencies change. The CI
# deploy-connect workflow regenerates manifest.json automatically.

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  message("Installing rsconnect ...")
  install.packages("rsconnect", quiet = TRUE)
}

## ---- manifest.json --------------------------------------------------------
# Honours .rscignore, so dev/, tests/, vignettes/ etc. are excluded.
rsconnect::writeManifest(appDir = ".", appPrimaryDoc = "app.R")
cat("Wrote manifest.json\n")

## ---- renv.lock (runtime scope) -------------------------------------------
# Scan only what the *deployed app* runs (entry point + package code + report
# template); the cellboard package itself is sourced from R/, not installed.
runtime_paths <- c("app.R", "R", "inst/report")
deps <- unique(renv::dependencies(runtime_paths, quiet = TRUE)$Package)
deps <- setdiff(deps, "cellboard")

renv::snapshot(
  project  = ".",
  lockfile = "renv.lock",
  packages = deps,
  prompt   = FALSE
)
cat("Wrote renv.lock with", length(deps), "top-level runtime packages:\n  ",
    paste(sort(deps), collapse = ", "), "\n")
