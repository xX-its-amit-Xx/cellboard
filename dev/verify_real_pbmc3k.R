#!/usr/bin/env Rscript
# End-to-end verification of the cellboard pipeline on the REAL public 10x
# pbmc3k dataset (downloaded by the caller to dev/cache/). Exercises every step
# a user / the cookbook performs and prints real numbers at each stage.

suppressPackageStartupMessages(library(ggplot2))
ok   <- function(...) cat(sprintf("[ OK ] %s\n", paste0(...)))
step <- function(...) cat(sprintf("\n=== %s ===\n", paste0(...)))
t0   <- function() proc.time()[["elapsed"]]

invisible(lapply(sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE)), source))

parent_dir <- "dev/cache/pbmc3k_extracted/filtered_gene_bc_matrices"
stopifnot(dir.exists(parent_dir))

step("1. Detect + load (pointing at the PARENT dir to test nested resolution)")
stopifnot(detect_input_type(parent_dir) == "10x_dir")
ok("detect_input_type() -> 10x_dir")
ti <- t0()
sce <- read_input(parent_dir)              # auto -> 10x, descends into hg19/
ok(sprintf("read_input(): %d genes x %d cells (%.1fs)",
           nrow(sce), ncol(sce), t0() - ti))
ok(sprintf("rowData cols: %s",
           paste(colnames(SummarizedExperiment::rowData(sce)), collapse = ", ")))
ok(sprintf("first rownames look like Ensembl IDs: %s",
           paste(head(rownames(sce), 2), collapse = ", ")))

step("2. Compute QC metrics on all ~2,700 cells (incl. scDblFinder)")
ti <- t0()
sce <- compute_qc(sce)
ok(sprintf("compute_qc(): %.1fs", t0() - ti))
cd <- as.data.frame(SummarizedExperiment::colData(sce))
meta <- S4Vectors::metadata(sce)$cellboard$qc
ok(sprintf("mito genes matched: %d | doublets called: %d",
           meta$n_mito_genes, sum(cd$doublet_class == "doublet", na.rm = TRUE)))
print(round(sapply(cd[, c("nUMI", "nGene", "pct_mito", "doublet_score")], summary), 2))

step("3. Threshold + filter (Seurat-tutorial-style cutoffs)")
th <- qc_thresholds(min_genes = 200, max_genes = 2500, max_mito = 5,
                    remove_doublet_class = TRUE)
res <- apply_qc_filter(sce, th)
ok(sprintf("retained %d / %d cells (%.1f%%), discarded %d",
           res$n_keep, res$n_total, 100 * res$n_keep / res$n_total, res$n_discard))
cat("Cells failing each filter:\n")
print(summarise_qc(res)$by_reason)

step("4. Export filtered SingleCellExperiment")
out_rds <- "dev/cache/pbmc3k_filtered_sce.rds"
saveRDS(res$filtered, out_rds)
chk <- readRDS(out_rds)
ok(sprintf("wrote %s: %d genes x %d cells (%.1f MB)",
           out_rds, nrow(chk), ncol(chk), file.info(out_rds)[["size"]] / 1048576))

step("5. Plots on real data")
dir.create("dev/cache", showWarnings = FALSE)
p <- qc_scatter(res$sce, "nUMI", "nGene", colour = "status", thresholds = th)
ggsave("dev/cache/pbmc3k_scatter.png", p, width = 7, height = 4.5, dpi = 110)
ok("wrote dev/cache/pbmc3k_scatter.png")

step("6. Render reproducible HTML report")
ti <- t0()
rep_path <- render_qc_report(res$sce, th, output_file = "pbmc3k_qc_report.html",
                             output_dir = "dev/cache", sample_name = "pbmc3k (real 10x)")
ok(sprintf("rendered %s (%.0f KB, %.1fs)",
           rep_path, file.info(rep_path)[["size"]] / 1024, t0() - ti))

step("7. Seurat interop round-trip")
suppressWarnings(suppressMessages({
  so  <- SeuratObject::CreateSeuratObject(
    counts = SummarizedExperiment::assay(sce, "counts"))
  sce2 <- as_sce(so)
  sce2 <- compute_qc(sce2, run_doublets = FALSE)
  keep <- apply_qc_filter(sce2, th)$keep
  so_qc <- so[, keep]
}))
ok(sprintf("Seurat -> SCE -> QC -> back to Seurat: %d -> %d cells",
           ncol(so), ncol(so_qc)))

cat("\nALL REAL-DATA STEPS PASSED\n")
