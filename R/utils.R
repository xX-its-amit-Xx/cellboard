# Internal helpers shared across the package. Everything here is fully
# namespace-qualified so the files can either be loaded as an installed package
# or sourced directly into the global environment by the Connect entry point
# (app.R) without behaviour changing.

# NULL-coalescing / empty-coalescing operator.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Class predicates ----------------------------------------------------------

is_sce <- function(x) methods::is(x, "SingleCellExperiment")

is_seurat <- function(x) methods::is(x, "Seurat")

# colData helpers -----------------------------------------------------------

# colData() as a plain data.frame, preserving cell names as row names.
coldata_df <- function(sce) {
  cd <- SummarizedExperiment::colData(sce)
  as.data.frame(cd, optional = TRUE, stringsAsFactors = FALSE)
}

# Add or overwrite a single colData column, preserving everything else. The
# replacement function returns a *new* object, so its result must be returned.
add_coldata <- function(sce, name, value) {
  cd <- SummarizedExperiment::colData(sce)
  cd[[name]] <- value
  SummarizedExperiment::`colData<-`(sce, value = cd)
}

# Attach a value under metadata()$cellboard, keeping other metadata intact.
set_cb_meta <- function(sce, key, value) {
  md <- S4Vectors::metadata(sce)
  if (is.null(md$cellboard)) md$cellboard <- list()
  md$cellboard[[key]] <- value
  S4Vectors::`metadata<-`(sce, value = md)
}

# Feature handling ----------------------------------------------------------

# Best guess at gene symbols, used for mitochondrial-gene pattern matching.
# Many 10x objects carry Ensembl IDs as row names and the symbol in rowData,
# so we look there first and fall back to the row names.
feature_symbols <- function(sce) {
  rd <- SummarizedExperiment::rowData(sce)
  cand <- c("Symbol", "symbol", "gene_symbol", "gene_symbols", "gene_name",
            "gene_names", "feature_name", "SYMBOL", "Gene", "gene", "name")
  hit <- intersect(cand, colnames(rd))
  sym <- if (length(hit)) as.character(rd[[hit[1L]]]) else rownames(sce)
  if (is.null(sym)) sym <- as.character(seq_len(nrow(sce)))
  # fall back to row names for any missing symbols
  rn <- rownames(sce) %||% as.character(seq_len(nrow(sce)))
  sym[is.na(sym) | sym == ""] <- rn[is.na(sym) | sym == ""]
  sym
}

# Validation ----------------------------------------------------------------

# Ensure the object is an SCE with a usable integer/numeric counts assay,
# coercing the sole assay to "counts" when it is unambiguous.
ensure_counts <- function(sce) {
  if (!is_sce(sce)) {
    stop("Expected a SingleCellExperiment, got ", class(sce)[1L], ".", call. = FALSE)
  }
  an <- SummarizedExperiment::assayNames(sce)
  if (!"counts" %in% an) {
    if (length(an) == 1L) {
      names(SummarizedExperiment::assays(sce))[1L] <- "counts"
    } else if (length(an) == 0L) {
      stop("The object has no assays; expected a 'counts' matrix.", call. = FALSE)
    } else {
      stop("No 'counts' assay found. Available assays: ",
           paste(an, collapse = ", "), ".", call. = FALSE)
    }
  }
  sce
}

# Locate a file shipped under inst/, working both when installed (system.file)
# and when the package tree is merely sourced (development / Connect bundle).
cb_system_file <- function(...) {
  p <- system.file(..., package = "cellboard")
  if (nzchar(p) && file.exists(p)) return(p)
  # development fallback: resolve relative to the package source tree
  local <- file.path("inst", ...)
  if (file.exists(local)) return(normalizePath(local, mustWork = FALSE))
  # second fallback: relative to this file's directory (app.R sourcing case)
  ""
}

# Compact human-readable integer, e.g. 12,345.
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
