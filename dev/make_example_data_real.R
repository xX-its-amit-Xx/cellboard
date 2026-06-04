#!/usr/bin/env Rscript
# Replace inst/extdata/ with a stratified 250-cell subset of the REAL
# 10x pbmc3k dataset. Stratified by mito% bin so the QC sliders have
# genuine structure: a healthy bulk, a stressed tail, and a few
# high-library-size outliers (likely doublets for scDblFinder to find).
#
# Run from the package root after dev/verify_real_pbmc3k.R has confirmed
# the pbmc3k download is present in dev/cache/.

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(DropletUtils)
})

set.seed(42)

pbmc_dir <- "dev/cache/pbmc3k_extracted/filtered_gene_bc_matrices/hg19"
out_dir  <- "inst/extdata"
tenx_out <- file.path(out_dir, "pbmc_small_10x")
rds_out  <- file.path(out_dir, "pbmc_small_sce.rds")

if (!dir.exists(pbmc_dir))
  stop("Run dev/verify_real_pbmc3k.R first to download the pbmc3k cache.")

cat("Reading full pbmc3k (32 738 genes x 2 700 cells)...\n")
full <- DropletUtils::read10xCounts(pbmc_dir, col.names = TRUE)
# Rename to standardise: DropletUtils puts ID + Symbol in rowData
if (!"Symbol" %in% colnames(SummarizedExperiment::rowData(full))) {
  rd <- SummarizedExperiment::rowData(full)
  colnames(rd)[1] <- "ID"
  colnames(rd)[2] <- "Symbol"
  SummarizedExperiment::rowData(full) <- rd
}
cnts <- SummarizedExperiment::assay(full, "counts")
symbols <- SummarizedExperiment::rowData(full)$Symbol
cat(sprintf("  genes: %d | cells: %d\n", nrow(full), ncol(full)))

# ---- Per-cell QC (direct, no Bioc deprecated path) -------------------------
mito_idx <- which(grepl("^MT-", symbols))
nUMI     <- Matrix::colSums(cnts)
nGene    <- Matrix::colSums(cnts > 0)
pct_mito <- 100 * Matrix::colSums(cnts[mito_idx, , drop = FALSE]) / nUMI
cat(sprintf("  mito genes: %d | median pct_mito: %.2f%%\n",
            length(mito_idx), median(pct_mito)))

# ---- Stratified sampling ---------------------------------------------------
# Bins: normal (mito < 3%), mild stress (3-8%), high mito (>8%),
#        high UMI outlier (top 3% nUMI, likely doublets)
high_umi   <- nUMI > quantile(nUMI, 0.97)
high_mito  <- pct_mito > 8  & !high_umi
mid_mito   <- pct_mito >= 3 & pct_mito <= 8 & !high_umi
normal     <- pct_mito < 3 & !high_umi

targets <- list(
  normal    = list(cells = which(normal),    n = 190),
  mid_mito  = list(cells = which(mid_mito),  n =  25),
  high_mito = list(cells = which(high_mito), n =  20),
  high_umi  = list(cells = which(high_umi),  n =  15)
)
cat(sprintf("  bin sizes -> normal:%d  mid_mito:%d  high_mito:%d  high_umi:%d\n",
            sum(normal), sum(mid_mito), sum(high_mito), sum(high_umi)))

selected <- unlist(lapply(targets, function(b) {
  sample(b$cells, min(b$n, length(b$cells)))
}))
selected <- sort(selected)
cat(sprintf("  selected: %d cells\n", length(selected)))

# ---- Subset cells + genes --------------------------------------------------
sub_cnts <- cnts[, selected, drop = FALSE]
# Keep genes detected in >= 2 of the selected cells; always retain all mito genes
detected_any <- Matrix::rowSums(sub_cnts > 0)
keep_genes   <- (detected_any >= 2) | seq_len(nrow(sub_cnts)) %in% mito_idx
sub_cnts <- sub_cnts[keep_genes, , drop = FALSE]
sub_ids  <- SummarizedExperiment::rowData(full)$ID[keep_genes]
sub_sym  <- symbols[keep_genes]
cat(sprintf("  -> %d genes x %d cells after gene filter\n",
            nrow(sub_cnts), ncol(sub_cnts)))

# ---- Build SCE -------------------------------------------------------------
sce <- SingleCellExperiment(
  assays  = list(counts = sub_cnts),
  rowData = S4Vectors::DataFrame(ID = sub_ids, Symbol = sub_sym),
  colData = S4Vectors::DataFrame(
    Barcode = colnames(sub_cnts),
    source  = "10x pbmc3k subset (real)"
  )
)
rownames(sce) <- sub_ids

# ---- Write -----------------------------------------------------------------
saveRDS(sce, rds_out)
if (dir.exists(tenx_out)) unlink(tenx_out, recursive = TRUE)
DropletUtils::write10xCounts(
  path        = tenx_out,
  x           = sub_cnts,
  barcodes    = colnames(sub_cnts),
  gene.id     = sub_ids,
  gene.symbol = sub_sym,
  gene.type   = "Gene Expression",
  version     = "3",
  overwrite   = TRUE
)
cat(sprintf("Wrote:\n  %s\n  %s/\n", rds_out, tenx_out))
cat(sprintf("  rds size: %.0f KB\n", file.info(rds_out)[["size"]] / 1024))
cat(sprintf("  10x files: %s\n",
            paste(list.files(tenx_out), collapse = ", ")))

# ---- Sanity check: source the package and run the pipeline -----------------
cat("\nRunning QC pipeline on the new example...\n")
invisible(lapply(sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE)), source))

qc <- compute_qc(sce, run_doublets = TRUE)
res <- apply_qc_filter(qc, default_thresholds(qc))
cat(sprintf("QC sanity check: %d/%d cells retained\n", res$n_keep, res$n_total))
cd <- as.data.frame(SummarizedExperiment::colData(qc))
cat("Median metrics:\n")
print(round(sapply(cd[, c("nUMI", "nGene", "pct_mito", "doublet_score")], median), 2))
