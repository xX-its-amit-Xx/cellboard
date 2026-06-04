suppressPackageStartupMessages({ library(cellboard); library(ggplot2) })

sce <- example_sce()
sce <- compute_qc(sce, run_doublets = TRUE)
res <- apply_qc_filter(sce, default_thresholds(sce))
cat(sprintf("Real example: %d genes x %d cells, %d/%d retained\n",
            nrow(sce), ncol(sce), res$n_keep, res$n_total))

p1 <- qc_scatter(res$sce, "nUMI", "nGene", colour = "status",
                 thresholds = default_thresholds(sce))
ggplot2::ggsave("man/figures/qc-scatter.png", p1, width = 7, height = 4.5, dpi = 110)

cd   <- as.data.frame(SummarizedExperiment::colData(sce))
long <- do.call(rbind, lapply(c("pct_mito", "nGene", "nUMI", "doublet_score"),
  function(m) data.frame(metric = m, value = cd[[m]],
                         status = ifelse(res$discard, "Discarded", "Retained"))))
long$metric <- factor(long$metric,
  levels = c("pct_mito", "nGene", "nUMI", "doublet_score"),
  labels = c("Mitochondrial %", "Genes / cell", "UMIs / cell", "Doublet score"))
p2 <- ggplot(long, aes("", value)) +
  geom_violin(fill = "grey92", colour = "grey55", scale = "width") +
  geom_jitter(aes(colour = status), width = 0.2, size = 0.7, alpha = 0.55) +
  scale_colour_manual(values = c(Retained = "#2C7FB8", Discarded = "#D95F02"), name = NULL) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL,
       title = "cellboard QC: real 10x pbmc3k subset (238 cells)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggplot2::ggsave("man/figures/qc-distributions.png", p2, width = 9, height = 3.6, dpi = 110)
cat("Updated man/figures/ with real pbmc3k data.\n")
