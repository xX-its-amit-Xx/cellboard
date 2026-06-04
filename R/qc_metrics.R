# The QC engine: compute per-cell metrics, define thresholds, and split cells
# into retained vs discarded. The filtering core (filter_cells_df) is a pure
# function over a data.frame so it can be unit-tested without Bioconductor.

#' Standardised QC metric column names
#'
#' The colData columns written by [compute_qc()] and read by the plots, the app
#' and the report.
#'
#' @return A named character vector mapping logical metric names to colData
#'   column names.
#' @export
qc_metric_cols <- function() .cb_cols

#' Compute per-cell QC metrics
#'
#' Adds standardised per-cell QC metrics to an SCE: UMIs (`nUMI`), genes
#' detected (`nGene`), mitochondrial percentage (`pct_mito`) and, optionally,
#' an [scDblFinder][scDblFinder::scDblFinder] doublet score (`doublet_score`)
#' and class (`doublet_class`).
#'
#' Mitochondrial genes are found by matching `mito_pattern` against gene symbols
#' (taken from a `Symbol`-like `rowData` column when present, otherwise the row
#' names), so it works whether your features are Ensembl IDs or symbols.
#'
#' @param sce A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'   with a `counts` assay.
#' @param mito_pattern Regular expression matching mitochondrial gene symbols.
#'   The default `"^[Mm][Tt]-"` matches human (`MT-`) and mouse (`mt-`).
#' @param run_doublets Logical; run scDblFinder (default `TRUE`). Falls back to
#'   `NA` doublet metrics if scDblFinder is unavailable or errors.
#' @param seed Random seed for the (stochastic) doublet step.
#' @param doublet_args Extra arguments passed to [scDblFinder::scDblFinder()].
#' @param verbose Emit a warning when no mitochondrial genes match.
#'
#' @return The input SCE with extra `colData` columns (see [qc_metric_cols()])
#'   and a `metadata(sce)$cellboard$qc` record.
#' @export
#'
#' @examples
#' sce <- example_sce()
#' sce <- compute_qc(sce, run_doublets = FALSE)
#' head(as.data.frame(SummarizedExperiment::colData(sce)))
compute_qc <- function(sce,
                       mito_pattern = "^[Mm][Tt]-",
                       run_doublets = TRUE,
                       seed = 100L,
                       doublet_args = list(),
                       verbose = TRUE) {
  sce <- ensure_counts(sce)

  syms <- feature_symbols(sce)
  mito <- which(grepl(mito_pattern, syms))
  n_mito <- length(mito)

  # Compute the per-cell metrics directly from the counts matrix. This is exact,
  # works on dense or sparse (Matrix) counts, and avoids depending on functions
  # that come and go across Bioconductor releases.
  cnts <- SummarizedExperiment::assay(sce, "counts")
  nUMI  <- Matrix::colSums(cnts)
  nGene <- Matrix::colSums(cnts > 0)
  mito_sum <- if (n_mito > 0L) Matrix::colSums(cnts[mito, , drop = FALSE]) else rep(0, ncol(sce))
  pct <- ifelse(nUMI > 0, 100 * mito_sum / nUMI, 0)
  if (n_mito == 0L && isTRUE(verbose)) {
    warning("No mitochondrial genes matched pattern '", mito_pattern,
            "'; pct_mito set to 0.", call. = FALSE)
  }

  sce <- add_coldata(sce, .cb_cols[["nUMI"]],  as.numeric(nUMI))
  sce <- add_coldata(sce, .cb_cols[["nGene"]], as.integer(nGene))
  sce <- add_coldata(sce, .cb_cols[["pct_mito"]], as.numeric(pct))

  # Doublet detection (stochastic; optional).
  ran <- FALSE
  if (isTRUE(run_doublets) && requireNamespace("scDblFinder", quietly = TRUE)) {
    res <- tryCatch({
      set.seed(seed)
      args <- c(list(sce = sce), doublet_args)
      suppressMessages(suppressWarnings(do.call(scDblFinder::scDblFinder, args)))
    }, error = function(e) {
      warning("scDblFinder failed (", conditionMessage(e),
              "); doublet metrics set to NA.", call. = FALSE)
      NULL
    })
    if (!is.null(res)) {
      dcd <- SummarizedExperiment::colData(res)
      ord <- match(colnames(sce), colnames(res))
      sce <- add_coldata(sce, .cb_cols[["doublet_score"]],
                         as.numeric(dcd$scDblFinder.score)[ord])
      sce <- add_coldata(sce, .cb_cols[["doublet_class"]],
                         as.character(dcd$scDblFinder.class)[ord])
      ran <- TRUE
    }
  }
  if (!ran) {
    sce <- add_coldata(sce, .cb_cols[["doublet_score"]], rep(NA_real_, ncol(sce)))
    sce <- add_coldata(sce, .cb_cols[["doublet_class"]], rep(NA_character_, ncol(sce)))
  }

  set_cb_meta(sce, "qc", list(
    mito_pattern = mito_pattern,
    n_mito_genes = n_mito,
    ran_doublets = ran,
    computed_at  = as.character(Sys.time())
  ))
}

#' Construct a set of QC thresholds
#'
#' Each threshold is either a single number (active) or `NULL` (inactive). A
#' cell is discarded if it fails *any* active threshold.
#'
#' @param min_genes,max_genes Minimum/maximum genes detected (`nGene`).
#' @param min_umi,max_umi Minimum/maximum UMIs (`nUMI`).
#' @param max_mito Maximum mitochondrial percentage (`pct_mito`).
#' @param max_doublet Maximum doublet score (`doublet_score`).
#' @param remove_doublet_class Logical; also discard cells whose
#'   `doublet_class` is `"doublet"`.
#'
#' @return An object of class `cb_thresholds` (a validated list).
#' @export
#'
#' @examples
#' qc_thresholds(min_genes = 200, max_mito = 10, remove_doublet_class = TRUE)
qc_thresholds <- function(min_genes = NULL, max_genes = NULL,
                          min_umi = NULL, max_umi = NULL,
                          max_mito = NULL, max_doublet = NULL,
                          remove_doublet_class = FALSE) {
  chk <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
      stop("`", nm, "` must be a single number or NULL.", call. = FALSE)
    }
    as.numeric(x)
  }
  th <- list(
    min_genes = chk(min_genes, "min_genes"),
    max_genes = chk(max_genes, "max_genes"),
    min_umi   = chk(min_umi, "min_umi"),
    max_umi   = chk(max_umi, "max_umi"),
    max_mito  = chk(max_mito, "max_mito"),
    max_doublet = chk(max_doublet, "max_doublet"),
    remove_doublet_class = isTRUE(remove_doublet_class)
  )
  if (!is.null(th$min_genes) && !is.null(th$max_genes) && th$min_genes > th$max_genes) {
    stop("min_genes (", th$min_genes, ") > max_genes (", th$max_genes, ").", call. = FALSE)
  }
  if (!is.null(th$min_umi) && !is.null(th$max_umi) && th$min_umi > th$max_umi) {
    stop("min_umi (", th$min_umi, ") > max_umi (", th$max_umi, ").", call. = FALSE)
  }
  structure(th, class = "cb_thresholds")
}

#' Sensible default QC thresholds
#'
#' Returns conservative starting thresholds. If an SCE with computed metrics is
#' supplied, an upper UMI bound is suggested from the data (median + 3 MADs on
#' the log scale) to catch likely multiplets.
#'
#' @param sce Optional SCE with QC metrics (from [compute_qc()]).
#'
#' @return A `cb_thresholds` object.
#' @export
default_thresholds <- function(sce = NULL) {
  th <- qc_thresholds(min_genes = 200, min_umi = 500, max_mito = 10,
                      remove_doublet_class = TRUE)
  if (!is.null(sce) && is_sce(sce)) {
    cd <- coldata_df(sce)
    if (.cb_cols[["nUMI"]] %in% names(cd)) {
      lu <- log10(cd[[.cb_cols[["nUMI"]]]] + 1)
      hi <- stats::median(lu, na.rm = TRUE) + 3 * stats::mad(lu, na.rm = TRUE)
      th$max_umi <- round(10^hi - 1)
    }
  }
  th
}

# Coerce list / cb_thresholds into a validated cb_thresholds.
as_thresholds <- function(x) {
  if (inherits(x, "cb_thresholds")) return(x)
  if (is.list(x)) return(do.call(qc_thresholds, x))
  stop("`thresholds` must be a cb_thresholds object or a named list.", call. = FALSE)
}

#' Decide which cells to keep, from a data.frame of metrics
#'
#' The pure filtering core. Operates on any data.frame that has the standardised
#' metric columns (see [qc_metric_cols()]); missing columns or `NA` values never
#' cause a cell to be discarded.
#'
#' @param df A data.frame with metric columns (`nUMI`, `nGene`, `pct_mito`,
#'   `doublet_score`, `doublet_class`).
#' @param thresholds A `cb_thresholds` object or a named list.
#'
#' @return A list with `keep` (logical), `discard` (logical) and `reasons`
#'   (a data.frame of per-filter logicals; `TRUE` = failed that filter).
#' @export
#'
#' @examples
#' df <- data.frame(nUMI = c(100, 5000), nGene = c(50, 1500), pct_mito = c(30, 5))
#' filter_cells_df(df, qc_thresholds(min_genes = 200, max_mito = 10))$keep
filter_cells_df <- function(df, thresholds = default_thresholds()) {
  th <- as_thresholds(thresholds)
  n <- nrow(df)
  col <- function(nm, default) if (nm %in% names(df)) df[[nm]] else rep(default, n)
  nGene  <- as.numeric(col(.cb_cols[["nGene"]], NA_real_))
  nUMI   <- as.numeric(col(.cb_cols[["nUMI"]], NA_real_))
  pct    <- as.numeric(col(.cb_cols[["pct_mito"]], NA_real_))
  dscore <- as.numeric(col(.cb_cols[["doublet_score"]], NA_real_))
  dclass <- as.character(col(.cb_cols[["doublet_class"]], NA_character_))

  # NA on a metric means "can't evaluate" -> not a failure.
  fail <- function(cond) { cond[is.na(cond)] <- FALSE; cond }
  reasons <- list()
  if (!is.null(th$min_genes))   reasons[["low_nGene"]]    <- fail(nGene < th$min_genes)
  if (!is.null(th$max_genes))   reasons[["high_nGene"]]   <- fail(nGene > th$max_genes)
  if (!is.null(th$min_umi))     reasons[["low_nUMI"]]     <- fail(nUMI  < th$min_umi)
  if (!is.null(th$max_umi))     reasons[["high_nUMI"]]    <- fail(nUMI  > th$max_umi)
  if (!is.null(th$max_mito))    reasons[["high_mito"]]    <- fail(pct   > th$max_mito)
  if (!is.null(th$max_doublet)) reasons[["doublet_score"]] <- fail(dscore > th$max_doublet)
  if (isTRUE(th$remove_doublet_class)) {
    reasons[["doublet_class"]] <- fail(dclass == "doublet")
  }

  reasons_df <- if (length(reasons)) {
    as.data.frame(reasons, stringsAsFactors = FALSE)
  } else {
    data.frame(row.names = seq_len(n))
  }
  discard <- if (ncol(reasons_df)) Reduce(`|`, reasons_df) else rep(FALSE, n)
  list(keep = !discard, discard = discard, reasons = reasons_df)
}

#' Apply QC thresholds to a SingleCellExperiment
#'
#' Annotates the SCE with `discard`/`qc_keep` columns and returns both the
#' annotated object and the filtered subset.
#'
#' @param sce An SCE with QC metrics (from [compute_qc()]).
#' @param thresholds A `cb_thresholds` object or named list.
#'
#' @return A list with elements `sce` (annotated), `filtered` (kept cells only),
#'   `keep`, `discard`, `reasons`, and the counts `n_total`, `n_keep`,
#'   `n_discard`.
#' @export
#'
#' @examples
#' sce <- compute_qc(example_sce(), run_doublets = FALSE)
#' res <- apply_qc_filter(sce, qc_thresholds(min_genes = 200, max_mito = 15))
#' res$n_keep
apply_qc_filter <- function(sce, thresholds = default_thresholds()) {
  sce <- ensure_counts(sce)
  th <- as_thresholds(thresholds)
  res <- filter_cells_df(coldata_df(sce), th)
  sce <- add_coldata(sce, "discard", res$discard)
  sce <- add_coldata(sce, "qc_keep", res$keep)
  sce <- set_cb_meta(sce, "thresholds", unclass(th))
  list(
    sce       = sce,
    filtered  = sce[, res$keep, drop = FALSE],
    keep      = res$keep,
    discard   = res$discard,
    reasons   = res$reasons,
    n_total   = ncol(sce),
    n_keep    = sum(res$keep),
    n_discard = sum(res$discard)
  )
}

#' Summarise a QC filtering result
#'
#' @param x An SCE (metrics computed), or the result of [apply_qc_filter()] /
#'   [filter_cells_df()].
#' @param thresholds Thresholds to use when `x` is an SCE.
#'
#' @return A list with `totals` (a data.frame: input/retained/discarded/percent)
#'   and `by_reason` (cells failing each active filter).
#' @export
summarise_qc <- function(x, thresholds = NULL) {
  res <- if (is_sce(x)) apply_qc_filter(x, thresholds %||% default_thresholds()) else x
  keep <- res$keep
  discard <- res$discard
  reasons <- res$reasons
  n_total <- length(keep)
  n_keep <- sum(keep)
  n_discard <- sum(discard)
  pct <- if (n_total) round(100 * n_keep / n_total, 1) else NA_real_

  lab <- c(low_nGene = "Too few genes", high_nGene = "Too many genes",
           low_nUMI = "Too few UMIs", high_nUMI = "Too many UMIs",
           high_mito = "High mitochondrial %", doublet_score = "High doublet score",
           doublet_class = "Called doublet (scDblFinder)")
  by_reason <- if (!is.null(reasons) && ncol(reasons)) {
    data.frame(
      filter        = names(reasons),
      description   = ifelse(is.na(lab[names(reasons)]), names(reasons), unname(lab[names(reasons)])),
      cells_failing = as.integer(vapply(reasons, sum, integer(1))),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    data.frame(filter = character(0), description = character(0),
               cells_failing = integer(0))
  }
  totals <- data.frame(
    metric = c("Input cells", "Retained", "Discarded", "Percent retained"),
    value  = c(n_total, n_keep, n_discard, pct),
    stringsAsFactors = FALSE
  )
  list(totals = totals, by_reason = by_reason)
}
