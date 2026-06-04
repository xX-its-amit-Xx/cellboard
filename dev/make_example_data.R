#!/usr/bin/env Rscript
# Generate the bundled example dataset: a small, deterministic, PBMC-style
# simulation with real marker + mitochondrial gene symbols and deliberately
# planted low-quality and doublet populations, so cellboard's QC sliders have
# clear structure to act on. Writes both a 10x v3 directory and an .rds SCE to
# inst/extdata/. This is SIMULATED data (no licensing constraints); the cookbook
# shows the identical workflow on the real 10x pbmc3k dataset.

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(DropletUtils)
})

set.seed(2024)

out_dir  <- "inst/extdata"
tenx_dir <- file.path(out_dir, "pbmc_small_10x")
rds_path <- file.path(out_dir, "pbmc_small_sce.rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Gene panel -----------------------------------------------------------
markers <- list(
  Tcell    = c("CD3D", "CD3E", "CD3G", "IL7R", "CCR7", "LDHB", "TCF7", "CD27"),
  CD8T     = c("CD8A", "CD8B", "GZMK", "CCL5", "NKG7"),
  Bcell    = c("CD79A", "CD79B", "MS4A1", "CD19", "HLA-DRA", "CD74"),
  Mono     = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1", "FCGR3A", "MS4A7"),
  NK       = c("GNLY", "NKG7", "KLRD1", "KLRF1", "NCAM1", "PRF1"),
  DC       = c("FCER1A", "CST3", "CLEC10A"),
  Platelet = c("PPBP", "PF4", "GP9", "ITGA2B")
)
housekeeping <- c("ACTB", "GAPDH", "B2M", "TMSB4X", "FTL", "FTH1", "TPT1", "MALAT1")
ribo <- c(paste0("RPS", c(2:6, 8, 11, 15, 18, 27)), paste0("RPL", c(3, 10, 13, 21, 23, 30, 34)))
mito <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
          "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB")

core_symbols <- unique(c(unlist(markers), housekeeping, ribo, mito))
n_genes <- 1200L
filler  <- sprintf("GENE%04d", seq_len(n_genes - length(core_symbols)))
symbols <- c(core_symbols, filler)
gene_ids <- sprintf("ENSG%011d", seq_len(n_genes))
names(symbols) <- symbols

# Per-gene heterogeneity (mean ~1) applied to the background rate only.
gene_var <- stats::rgamma(n_genes, shape = 2, rate = 2)
names(gene_var) <- symbols

DEPTH <- 10  # global multiplier -> realistic library sizes

# Baseline per-gene expression rate (lambda) for a cell of a given type.
base_lambda <- function(type) {
  lam <- rep(0.03, n_genes) * gene_var
  names(lam) <- symbols
  lam[housekeeping] <- 6
  lam[ribo] <- 4
  lam[mito] <- 0.8
  for (m in markers[[type]]) if (m %in% names(lam)) lam[m] <- 12
  other <- sample(setdiff(names(markers), type), 1)            # leaky lineage
  for (m in markers[[other]]) if (m %in% names(lam)) lam[m] <- 1.5
  lam
}

## ---- Cell design ----------------------------------------------------------
n_good    <- 270L  # healthy
n_lowq    <- 30L   # low-quality: high mito %, low complexity
n_doublet <- 20L   # planted doublets: ~2x library, mixed lineage
n_cells   <- n_good + n_lowq + n_doublet
cell_types <- names(markers)

type_good <- sample(cell_types, n_good, replace = TRUE,
                    prob = c(0.30, 0.12, 0.15, 0.20, 0.12, 0.05, 0.06))
type_lowq <- sample(cell_types, n_lowq, replace = TRUE)

counts <- matrix(0L, nrow = n_genes, ncol = n_cells,
                 dimnames = list(gene_ids, sprintf("cell_%03d", seq_len(n_cells))))

# Healthy cells
size_good <- stats::rlnorm(n_good, 0, 0.30)
for (j in seq_len(n_good)) {
  counts[, j] <- stats::rpois(n_genes, base_lambda(type_good[j]) * DEPTH * size_good[j])
}
# Low-quality cells: suppress non-mito, inflate mito -> high pct_mito, low nGene
size_lowq <- stats::rlnorm(n_lowq, log(0.6), 0.25)
for (k in seq_len(n_lowq)) {
  j <- n_good + k
  lam <- base_lambda(type_lowq[k]) * DEPTH * size_lowq[k]
  is_mito <- names(lam) %in% mito
  lam[!is_mito] <- lam[!is_mito] * 0.20
  lam[is_mito]  <- 0.8 * DEPTH * stats::runif(1, 4, 9)   # strong mitochondrial signal
  counts[, j] <- stats::rpois(n_genes, lam)
}
# Doublets: sum of two random lineages, larger library
size_db <- stats::rlnorm(n_doublet, log(1.5), 0.20)
for (d in seq_len(n_doublet)) {
  j <- n_good + n_lowq + d
  pair <- sample(cell_types, 2)
  lam <- (base_lambda(pair[1]) + base_lambda(pair[2])) * DEPTH * size_db[d]
  counts[, j] <- stats::rpois(n_genes, lam)
}

# Keep genes seen in at least one cell, but always retain mito genes.
keep_genes <- (Matrix::rowSums(counts) > 0) | (symbols %in% mito)
counts   <- counts[keep_genes, , drop = FALSE]
gene_ids <- gene_ids[keep_genes]
symbols  <- symbols[keep_genes]
sparse_counts <- as(counts, "CsparseMatrix")

## ---- Assemble SCE ---------------------------------------------------------
sce <- SingleCellExperiment(
  assays  = list(counts = sparse_counts),
  rowData = DataFrame(ID = gene_ids, Symbol = unname(symbols)),
  colData = DataFrame(
    Barcode = colnames(sparse_counts),
    truth   = rep(c("cell", "low_quality", "doublet"),
                  c(n_good, n_lowq, n_doublet))   # hidden ground truth
  )
)
rownames(sce) <- gene_ids
cat(sprintf("Simulated SCE: %d genes x %d cells\n", nrow(sce), ncol(sce)))

## ---- Write outputs --------------------------------------------------------
saveRDS(sce, rds_path)
if (dir.exists(tenx_dir)) unlink(tenx_dir, recursive = TRUE)
DropletUtils::write10xCounts(
  path = tenx_dir, x = sparse_counts, barcodes = colnames(sparse_counts),
  gene.id = gene_ids, gene.symbol = unname(symbols),
  gene.type = "Gene Expression", version = "3", overwrite = TRUE
)
cat("Wrote:\n  ", rds_path, "\n  ", tenx_dir, "/ (10x v3)\n", sep = "")
cat("Files in 10x dir: ", paste(list.files(tenx_dir), collapse = ", "), "\n")

## ---- Sanity check: run the QC pipeline ------------------------------------
invisible(lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))
qc  <- compute_qc(sce, run_doublets = TRUE)
res <- apply_qc_filter(qc, default_thresholds())
cat(sprintf("\nQC sanity check: %d/%d cells retained under default thresholds\n",
            res$n_keep, res$n_total))
cd <- as.data.frame(SummarizedExperiment::colData(qc))
cat("Median metrics by ground-truth group:\n")
print(aggregate(cbind(nUMI, nGene, pct_mito, doublet_score) ~ truth, data = cd, FUN = median))
cat("\nDiscarded vs ground truth:\n")
print(table(truth = cd$truth, discarded = res$discard))
