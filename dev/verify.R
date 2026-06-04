#!/usr/bin/env Rscript
# End-to-end smoke test of the cellboard sources (no install required):
#   sources R/, runs the QC pipeline, drives the Shiny server with testServer(),
#   renders the report, and writes README figures to man/figures/.

suppressPackageStartupMessages({
  library(ggplot2)
})

ok <- function(msg) cat(sprintf("[ OK ] %s\n", msg))
step <- function(msg) cat(sprintf("\n=== %s ===\n", msg))

step("Source R/")
invisible(lapply(sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE)), source))
ok("sourced all R/ files")

step("I/O")
sce <- example_sce()
stopifnot(is_sce(sce), "counts" %in% SummarizedExperiment::assayNames(sce))
ok(sprintf("example_sce(): %d genes x %d cells", nrow(sce), ncol(sce)))
sce10x <- read_10x(example_sce(as = "10x_path"))
stopifnot(ncol(sce10x) == ncol(sce))
ok("read_10x() round-trips the bundled 10x directory")
stopifnot(detect_input_type(example_sce(as = "10x_path")) == "10x_dir")
ok("detect_input_type() identifies the 10x directory")

step("QC engine")
qc <- compute_qc(sce, run_doublets = TRUE)
cd <- as.data.frame(SummarizedExperiment::colData(qc))
stopifnot(all(c("nUMI", "nGene", "pct_mito", "doublet_score") %in% names(cd)))
ok(sprintf("compute_qc(): median mito %.1f%%, doublets flagged %d",
           median(cd$pct_mito), sum(cd$doublet_class == "doublet", na.rm = TRUE)))

res <- apply_qc_filter(qc, default_thresholds(qc))
ok(sprintf("apply_qc_filter(): %d/%d retained", res$n_keep, res$n_total))
stopifnot(res$n_keep < res$n_total, res$n_keep > 0)

# stricter thresholds must retain no more cells (monotonicity)
strict <- apply_qc_filter(qc, qc_thresholds(min_genes = 500, max_mito = 5))
stopifnot(strict$n_keep <= res$n_keep)
ok("stricter thresholds retain fewer cells (monotone)")

step("Pure filter core (no Bioc)")
df <- data.frame(nUMI = c(100, 5000, 8000), nGene = c(40, 1500, 2200),
                 pct_mito = c(60, 4, 5),
                 doublet_class = c("singlet", "singlet", "doublet"))
fk <- filter_cells_df(df, qc_thresholds(min_genes = 200, max_mito = 10,
                                        remove_doublet_class = TRUE))$keep
stopifnot(identical(fk, c(FALSE, TRUE, FALSE)))
ok("filter_cells_df() drops low-quality and doublet rows as expected")

step("Plots")
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
p_scatter <- qc_scatter(res$sce, "nUMI", "nGene", colour = "status",
                        thresholds = default_thresholds(qc))
ggsave("man/figures/qc-scatter.png", p_scatter, width = 7, height = 4.5, dpi = 110)
ok("wrote man/figures/qc-scatter.png")

# faceted distributions figure
long <- do.call(rbind, lapply(c("pct_mito", "nGene", "nUMI", "doublet_score"),
  function(m) data.frame(metric = m, value = cd[[m]],
                         status = ifelse(res$discard, "Discarded", "Retained"))))
long$metric <- factor(long$metric,
  levels = c("pct_mito", "nGene", "nUMI", "doublet_score"),
  labels = c("Mitochondrial %", "Genes / cell", "UMIs / cell", "Doublet score"))
p_dist <- ggplot(long, aes("", value)) +
  geom_violin(fill = "grey92", colour = "grey55", scale = "width") +
  geom_jitter(aes(colour = status), width = 0.2, size = 0.5, alpha = 0.5) +
  scale_colour_manual(values = c(Retained = "#2C7FB8", Discarded = "#D95F02"),
                      name = NULL) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, title = "cellboard QC metric distributions") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("man/figures/qc-distributions.png", p_dist, width = 9, height = 3.6, dpi = 110)
ok("wrote man/figures/qc-distributions.png")

step("Report render")
out <- render_qc_report(res$sce, default_thresholds(qc),
                        output_file = "verify_report.html",
                        output_dir = "dev", quiet = TRUE)
stopifnot(file.exists(out), file.info(out)$size > 5000)
ok(sprintf("rendered report: %s (%.0f KB)", out, file.info(out)$size / 1024))

step("Shiny server (testServer)")
shiny::testServer(app_server, {
  session$setInputs(load_example = 1)
  session$setInputs(mito_pattern = "^[Mm][Tt]-", run_doublets = TRUE)
  session$setInputs(compute = 1)
  session$setInputs(max_mito = 10)
  out_keep <- filt()$n_keep
  out_total <- filt()$n_total
  cat(sprintf("  server filt(): %d/%d retained\n", out_keep, out_total))
  stopifnot(out_keep > 0, out_keep < out_total)
})
ok("testServer drove load -> compute -> filter")

step("UI builds")
ui <- app_ui()
stopifnot(inherits(ui, "shiny.tag.list") || inherits(ui, "shiny.tag") ||
          inherits(ui, "bslib_page"))
ok("app_ui() constructs")

cat("\nALL VERIFICATION STEPS PASSED\n")
