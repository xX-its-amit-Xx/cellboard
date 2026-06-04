# Input/output: turn whatever the user has (10x, SingleCellExperiment, Seurat)
# into a validated SingleCellExperiment with a "counts" assay.

#' Detect the type of a single-cell input path
#'
#' Inspects a file or directory and guesses how it should be read.
#'
#' @param path Path to a file or directory.
#'
#' @return One of `"10x_dir"`, `"10x_h5"`, `"rds"`, or throws an error if the
#'   path does not look like a supported input.
#' @export
#'
#' @examples
#' detect_input_type(cellboard:::cb_system_file("extdata", "pbmc_small_10x"))
detect_input_type <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a single non-empty string.", call. = FALSE)
  }
  if (dir.exists(path)) {
    files <- list.files(path)
    has_mtx <- any(grepl("^matrix\\.mtx(\\.gz)?$", files))
    has_bc  <- any(grepl("^barcodes\\.tsv(\\.gz)?$", files))
    has_feat <- any(grepl("^(features|genes)\\.tsv(\\.gz)?$", files))
    if (has_mtx && has_bc && has_feat) return("10x_dir")
    # Cell Ranger nests these under filtered_feature_bc_matrix/, so recurse one level.
    sub <- list.dirs(path, recursive = FALSE)
    for (d in sub) {
      f <- list.files(d)
      if (any(grepl("^matrix\\.mtx(\\.gz)?$", f)) &&
          any(grepl("^barcodes\\.tsv(\\.gz)?$", f))) {
        return("10x_dir")
      }
    }
    stop("Directory '", path, "' does not contain a 10x matrix ",
         "(matrix.mtx[.gz], barcodes.tsv[.gz], features/genes.tsv[.gz]).",
         call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Path does not exist: ", path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  switch(ext,
    h5     = "10x_h5",
    hdf5   = "10x_h5",
    rds    = "rds",
    # a bare matrix.mtx[.gz] -> treat its directory as a 10x folder
    mtx    = "10x_dir",
    gz     = if (grepl("matrix\\.mtx\\.gz$", path)) "10x_dir" else
               stop("Unsupported gzip input: ", path, call. = FALSE),
    stop("Unsupported input extension '.", ext, "'. Supported: directory of ",
         "10x files, .h5/.hdf5 (10x HDF5), or .rds (SCE/Seurat).", call. = FALSE)
  )
}

#' Read a 10x Genomics matrix into a SingleCellExperiment
#'
#' Reads either a Cell Ranger output directory (`matrix.mtx[.gz]` +
#' `barcodes.tsv[.gz]` + `features.tsv[.gz]`) or a 10x HDF5 (`.h5`) file.
#'
#' @param path Path to the 10x directory or `.h5` file.
#' @param ... Passed to [DropletUtils::read10xCounts()].
#'
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment].
#' @export
read_10x <- function(path, ...) {
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    stop("Reading 10x input requires the 'DropletUtils' package.", call. = FALSE)
  }
  # If a bare matrix.mtx[.gz] was given, point at its directory.
  if (!dir.exists(path) && grepl("matrix\\.mtx(\\.gz)?$", path)) {
    path <- dirname(path)
  }
  # Cell Ranger often nests the matrix one or two levels down
  # (e.g. filtered_gene_bc_matrices/hg19/); descend to the real matrix dir.
  if (dir.exists(path)) path <- resolve_10x_dir(path)
  sce <- DropletUtils::read10xCounts(path, col.names = TRUE, ...)
  ensure_counts(sce)
}

# Find the directory that actually holds matrix.mtx[.gz], descending into
# subdirectories if necessary. Returns `path` unchanged if none is found
# (so read10xCounts can raise its own informative error).
resolve_10x_dir <- function(path) {
  has_matrix <- function(d) any(grepl("^matrix\\.mtx(\\.gz)?$", list.files(d)))
  if (has_matrix(path)) return(path)
  subs <- list.dirs(path, recursive = TRUE, full.names = TRUE)
  hit <- Filter(has_matrix, subs)
  if (length(hit)) hit[[1L]] else path
}

#' Read a SingleCellExperiment from an .rds file
#'
#' @param path Path to an `.rds` file containing a SingleCellExperiment (or an
#'   object coercible to one, e.g. a SummarizedExperiment or a count matrix).
#'
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment].
#' @export
read_sce <- function(path) {
  obj <- readRDS(path)
  as_sce(obj)
}

#' Read a Seurat object from an .rds file and convert it to a SingleCellExperiment
#'
#' @param path Path to an `.rds` file containing a Seurat object.
#' @param assay Seurat assay to pull counts from (default `"RNA"`); ignored if
#'   absent.
#'
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment].
#' @export
read_seurat <- function(path, assay = "RNA") {
  obj <- readRDS(path)
  if (!is_seurat(obj)) {
    stop("File '", path, "' does not contain a Seurat object (got ",
         class(obj)[1L], "). Use read_sce() instead.", call. = FALSE)
  }
  as_sce(obj, assay = assay)
}

#' Coerce a supported object to a SingleCellExperiment
#'
#' Accepts a SingleCellExperiment (returned unchanged), a Seurat object
#' (converted via [Seurat::as.SingleCellExperiment()]), a SummarizedExperiment,
#' or a (sparse) count matrix.
#'
#' @param x Object to coerce.
#' @param assay For Seurat input, the assay to convert (default `"RNA"`).
#'
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'   with a `counts` assay.
#' @export
as_sce <- function(x, assay = "RNA") {
  if (is_sce(x)) {
    return(ensure_counts(x))
  }
  if (is_seurat(x)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Converting Seurat input requires the 'Seurat' package.", call. = FALSE)
    }
    assays_present <- tryCatch(SeuratObject::Assays(x), error = function(e) NULL)
    use_assay <- if (!is.null(assays_present) && assay %in% assays_present) assay else NULL
    sce <- if (is.null(use_assay)) {
      Seurat::as.SingleCellExperiment(x)
    } else {
      Seurat::as.SingleCellExperiment(x, assay = use_assay)
    }
    return(ensure_counts(sce))
  }
  if (methods::is(x, "SummarizedExperiment")) {
    return(ensure_counts(methods::as(x, "SingleCellExperiment")))
  }
  if (methods::is(x, "Matrix") || is.matrix(x)) {
    m <- methods::as(x, "CsparseMatrix")
    if (is.null(colnames(m))) colnames(m) <- paste0("cell", seq_len(ncol(m)))
    if (is.null(rownames(m))) rownames(m) <- paste0("gene", seq_len(nrow(m)))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
    return(sce)
  }
  stop("Don't know how to coerce an object of class '", class(x)[1L],
       "' to a SingleCellExperiment.", call. = FALSE)
}

#' Load any supported single-cell input
#'
#' The high-level entry point used by the dashboard. Detects the input type
#' (or uses the one you supply) and returns a validated SingleCellExperiment.
#'
#' @param path Path to a 10x directory, a 10x `.h5` file, or an `.rds` holding a
#'   SingleCellExperiment or Seurat object.
#' @param type One of `"auto"` (default), `"10x"`, `"sce"`, or `"seurat"`.
#' @param ... Passed to the underlying reader.
#'
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'   with a `counts` assay.
#' @export
#'
#' @examples
#' \dontrun{
#'   sce <- read_input("filtered_feature_bc_matrix/")
#'   sce <- read_input("pbmc.rds", type = "seurat")
#' }
read_input <- function(path, type = c("auto", "10x", "sce", "seurat"), ...) {
  type <- match.arg(type)
  if (type == "auto") {
    detected <- detect_input_type(path)
    type <- switch(detected,
      "10x_dir" = "10x",
      "10x_h5"  = "10x",
      "rds"     = "rds"
    )
  }
  sce <- switch(type,
    "10x"    = read_10x(path, ...),
    "sce"    = read_sce(path),
    "seurat" = read_seurat(path, ...),
    "rds"    = read_sce(path),   # readRDS + class-based coercion
    stop("Unknown input type: ", type, call. = FALSE)
  )
  ensure_counts(sce)
}

#' Load the bundled example PBMC dataset
#'
#' Returns a stratified 238-cell subset of the **real** 10x Genomics pbmc3k
#' dataset (hg19, Cell Ranger v2 output), with 9,936 expressed genes. The
#' subset preserves the biologically realistic QC structure of the full dataset:
#' a healthy bulk (~2% mitochondrial reads), a stressed tail, and high-UMI
#' outliers, making it suitable for demonstrating cellboard's QC workflow.
#' Generated from the full 2,700-cell dataset by `dev/make_example_data_real.R`.
#'
#' @param as A character scalar: `"sce"` (default) returns a
#'   SingleCellExperiment; `"10x_path"` returns the path to the bundled 10x
#'   directory (useful for demonstrating [read_10x()]).
#'
#' @return A SingleCellExperiment, or a path string.
#' @export
#'
#' @examples
#' sce <- example_sce()
#' dim(sce)
example_sce <- function(as = c("sce", "10x_path")) {
  as <- match.arg(as)
  if (as == "10x_path") {
    p <- cb_system_file("extdata", "pbmc_small_10x")
    if (!nzchar(p)) stop("Bundled example 10x directory not found.", call. = FALSE)
    return(p)
  }
  p <- cb_system_file("extdata", "pbmc_small_sce.rds")
  if (nzchar(p)) return(ensure_counts(readRDS(p)))
  # Fall back to reading the bundled 10x directory if the rds is unavailable.
  p10 <- cb_system_file("extdata", "pbmc_small_10x")
  if (nzchar(p10)) return(read_10x(p10))
  stop("Bundled example dataset not found.", call. = FALSE)
}
