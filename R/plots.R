# Plot helpers. Each returns a ggplot object built from a small tidy data.frame
# (so column names are literal and no tidy-eval pronoun is needed). The Shiny
# app wraps these with plotly::ggplotly() for interactivity; the report uses
# them as static figures.

# Discrete colours for retained/discarded status.
.cb_status_cols <- c(Retained = "#2C7FB8", Discarded = "#D95F02", Cells = "#636363")

# Per-cell status from a colData data.frame.
status_vec <- function(cd) {
  if ("discard" %in% names(cd)) {
    ifelse(isTRUE_vec(cd$discard), "Discarded", "Retained")
  } else {
    rep("Cells", nrow(cd))
  }
}
isTRUE_vec <- function(x) !is.na(x) & x

# Active threshold values relevant to a metric, for drawing reference lines.
metric_threshold_lines <- function(metric, thresholds) {
  if (is.null(thresholds)) return(numeric(0))
  th <- as_thresholds(thresholds)
  v <- switch(metric,
    nGene         = c(th$min_genes, th$max_genes),
    nUMI          = c(th$min_umi, th$max_umi),
    pct_mito      = th$max_mito,
    doublet_score = th$max_doublet,
    numeric(0)
  )
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  # a 0-valued lower bound is no filter; also it breaks log axes -> drop it
  if (metric %in% c("nUMI", "nGene")) v <- v[v > 0]
  v
}

# Safe metric label lookup: returns the human label or the key itself.
lbl <- function(k) {
  if (!is.null(k) && length(k) == 1L && k %in% names(.cb_labels)) {
    unname(.cb_labels[[k]])
  } else {
    k
  }
}

# A shared minimal theme.
cb_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# A placeholder plot used when a metric is unavailable.
blank_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 5, colour = "grey40") +
    ggplot2::theme_void()
}

#' Violin plot of a QC metric
#'
#' @param sce An SCE with QC metrics (and optionally a `discard` column).
#' @param metric One of `"pct_mito"`, `"nGene"`, `"nUMI"`, `"doublet_score"`.
#' @param thresholds Optional thresholds to draw as dashed reference lines.
#' @param points Overlay jittered per-cell points (default `TRUE`).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' sce <- compute_qc(example_sce(), run_doublets = FALSE)
#' qc_violin(sce, "pct_mito", qc_thresholds(max_mito = 10))
qc_violin <- function(sce,
                      metric = c("pct_mito", "nGene", "nUMI", "doublet_score"),
                      thresholds = NULL,
                      points = TRUE) {
  metric <- match.arg(metric)
  cd <- coldata_df(sce)
  col <- .cb_cols[[metric]]
  if (!col %in% names(cd) || all(is.na(cd[[col]]))) {
    return(blank_plot(paste0(.cb_labels[[metric]], " not available")))
  }
  pd <- data.frame(
    value  = as.numeric(cd[[col]]),
    group  = .cb_labels[[metric]],
    status = status_vec(cd),
    stringsAsFactors = FALSE
  )
  pd <- pd[!is.na(pd$value), , drop = FALSE]
  # values <= 0 cannot be shown on a log axis; drop them to avoid -Inf warnings
  if (metric %in% c("nUMI", "nGene")) pd <- pd[pd$value > 0, , drop = FALSE]

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = group, y = value)) +
    ggplot2::geom_violin(fill = "grey92", colour = "grey55", scale = "width", width = 0.85)
  if (isTRUE(points)) {
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(colour = status),
      width = 0.18, height = 0, size = 0.7, alpha = 0.55
    ) +
      ggplot2::scale_colour_manual(values = .cb_status_cols, name = NULL)
  }
  for (ln in metric_threshold_lines(metric, thresholds)) {
    p <- p + ggplot2::geom_hline(yintercept = ln, linetype = "dashed",
                                 colour = "#B2182B", linewidth = 0.6)
  }
  if (metric %in% c("nUMI", "nGene")) {
    p <- p + ggplot2::scale_y_log10(labels = scales_comma())
  }
  p + ggplot2::labs(x = NULL, y = .cb_labels[[metric]],
                    title = .cb_labels[[metric]]) +
    cb_theme()
}

# Light wrapper so we don't hard-depend on the 'scales' package being attached.
scales_comma <- function() {
  function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' Scatter plot of two QC metrics
#'
#' Typically UMIs vs genes, coloured by retained/discarded status or by a
#' continuous metric. Threshold reference lines are drawn for the plotted axes.
#'
#' @param sce An SCE with QC metrics.
#' @param x,y Metric column names for the axes (default `"nUMI"`, `"nGene"`).
#' @param colour How to colour points: `"status"` (retained vs discarded),
#'   `"pct_mito"`, or `"doublet_score"`. Defaults to `"status"` when a `discard`
#'   column is present, otherwise `"pct_mito"`.
#' @param thresholds Optional thresholds to draw as reference lines.
#'
#' @return A ggplot object.
#' @export
qc_scatter <- function(sce,
                       x = "nUMI", y = "nGene",
                       colour = NULL,
                       thresholds = NULL) {
  cd <- coldata_df(sce)
  if (!all(c(x, y) %in% names(cd))) {
    return(blank_plot("Metrics not available"))
  }
  colour <- colour %||% (if ("discard" %in% names(cd)) "status" else "pct_mito")
  pd <- data.frame(x = as.numeric(cd[[x]]), y = as.numeric(cd[[y]]),
                   stringsAsFactors = FALSE)
  discrete <- identical(colour, "status")
  pd$col <- if (discrete) status_vec(cd) else as.numeric(cd[[colour]])
  # drop non-positive values on axes that will be log-scaled
  if (x %in% c("nUMI", "nGene")) pd <- pd[pd$x > 0, , drop = FALSE]
  if (y %in% c("nUMI", "nGene")) pd <- pd[pd$y > 0, , drop = FALSE]

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = y, colour = col)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.7)
  if (discrete) {
    p <- p + ggplot2::scale_colour_manual(values = .cb_status_cols, name = NULL)
  } else {
    p <- p + ggplot2::scale_colour_viridis_c(name = lbl(colour))
  }
  # vertical lines for x thresholds, horizontal for y thresholds
  for (ln in metric_threshold_lines(x, thresholds)) {
    p <- p + ggplot2::geom_vline(xintercept = ln, linetype = "dashed",
                                 colour = "#B2182B", linewidth = 0.5)
  }
  for (ln in metric_threshold_lines(y, thresholds)) {
    p <- p + ggplot2::geom_hline(yintercept = ln, linetype = "dashed",
                                 colour = "#B2182B", linewidth = 0.5)
  }
  if (x %in% c("nUMI", "nGene")) p <- p + ggplot2::scale_x_log10(labels = scales_comma())
  if (y %in% c("nUMI", "nGene")) p <- p + ggplot2::scale_y_log10(labels = scales_comma())
  p + ggplot2::labs(x = lbl(x), y = lbl(y),
                    title = paste0(lbl(y), " vs ", lbl(x))) +
    cb_theme()
}

#' Barcode-rank (knee) plot
#'
#' Cells ranked by total UMIs, on log-log axes; useful for eyeballing an
#' empty-droplet / low-count cutoff.
#'
#' @param sce An SCE with QC metrics.
#' @param thresholds Optional thresholds (draws `min_umi` as a reference line).
#'
#' @return A ggplot object.
#' @export
qc_knee <- function(sce, thresholds = NULL) {
  cd <- coldata_df(sce)
  col <- .cb_cols[["nUMI"]]
  if (!col %in% names(cd)) return(blank_plot("UMI counts not available"))
  totals <- sort(as.numeric(cd[[col]]), decreasing = TRUE)
  totals <- totals[totals > 0]
  pd <- data.frame(rank = seq_along(totals), umi = totals)
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = rank, y = umi)) +
    ggplot2::geom_line(colour = "#2C7FB8", linewidth = 0.9) +
    ggplot2::scale_x_log10(labels = scales_comma()) +
    ggplot2::scale_y_log10(labels = scales_comma())
  for (ln in metric_threshold_lines("nUMI", thresholds)) {
    p <- p + ggplot2::geom_hline(yintercept = ln, linetype = "dashed",
                                 colour = "#B2182B", linewidth = 0.6)
  }
  p + ggplot2::labs(x = "Barcode rank", y = "Total UMIs",
                    title = "Barcode-rank (knee) plot") +
    cb_theme()
}
