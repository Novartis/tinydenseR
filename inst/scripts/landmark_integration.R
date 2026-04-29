###############################################################################
# landmark_integration.R
#
# Systematic comparison: landmark-level integration (tinydenseR / RunTDR)
# vs cell-level integration (Harmony) on the Luecken et al. Immune_ALL_human
# benchmarking dataset.
#
# Dataset: https://figshare.com/ndownloader/files/25717328
#          Luecken et al., Nat Methods 2021
#          Immune_ALL_human.h5ad  (33,506 cells × 12,303 genes)
#          10 batches, 5 studies, 16 cell types (final_annotation)
#
# Copyright 2026 Novartis Biomedical Research Inc.
# SPDX-License-Identifier: MIT
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 0. Environment & reproducibility
# ═══════════════════════════════════════════════════════════════════════════════

set.seed(42L)

library(tinydenseR)
library(sparseMatrixStats)
library(BPCells)
library(rhdf5)          # h5ad reading (used by tinydenseR internals)
library(harmony)
library(symphony)
library(irlba)
library(Matrix)
library(RANN)           # exact kNN
library(igraph)         # graph connectivity metric
library(cluster)        # silhouette
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(patchwork)

# set wd to local dir
script.path <-
  grep(pattern = "^--file=", 
       x = commandArgs(), 
       value = TRUE,
       fixed = FALSE) |>
  gsub(pattern = "^--file=",
       replacement = "", 
       fixed = FALSE) |>
  (\(x)
   if(length(x = x) == 0) rstudioapi::getSourceEditorContext()$path else x
  )()

script.path |>
  dirname() |>
  setwd()

# ── Paths (explicit, no hidden state) ─────────────────────────────────────────
# Required data: Immune_ALL_human.h5ad (33,506 cells × 12,303 genes)
# Download from: https://figshare.com/ndownloader/files/25717328
# Place in:      derived_data/Immune_ALL_human.h5ad
#
# Example download command (R):
#   dir.create("derived_data", showWarnings = FALSE)
#   curl::curl_download(
#     "https://figshare.com/ndownloader/files/25717328",
#     file.path("derived_data", "Immune_ALL_human.h5ad"))

H5AD_PATH    <- file.path("derived_data", "Immune_ALL_human.h5ad")
DERIVED_DIR  <- "derived_data"
RESULTS_DIR  <- "results"
BPCELLS_DIR  <- file.path(DERIVED_DIR, "Immune_ALL_human_bpcells")

if (!dir.exists(DERIVED_DIR)) dir.create(DERIVED_DIR, recursive = TRUE)
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ── Global settings ───────────────────────────────────────────────────────────

SEED                <- 42L
N_HVG               <- 5000L      # highly variable genes
N_PC                <- 30L         # principal components
BATCH_VAR           <- "batch"     # column used for Harmony integration
CELLTYPE_VAR        <- "final_annotation"
PSEUDO_SAMPLE_VAR   <- "pseudo_sample"  # synthetic sample grouping (see §1)
N_PSEUDO_PER_BATCH  <- 3L          # random sub-samples per batch
PROP_LANDMARKS      <- 0.1         # 10% landmark sampling rate
N_THREADS           <- 1L

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 1. Data preparation (shared between branches)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Open the BPCells on-disk matrix, read cell metadata, and create pseudo-
# samples.
#
# WHY PSEUDO-SAMPLES:
#   tinydenseR requires a `.sample.var` column that defines the grouping unit
#   for (a) proportional landmark selection and (b) per-sample Symphony
#   projection. The Luecken Immune_ALL_human dataset has no natural sub-batch
#   sample structure: `sample_ID` (5 values) is 1:1 with `study` and maps
#   identically to the batch grouping. Using `.sample.var = "batch"` would
#   conflate the *sample grouping* (which controls landmark selection and
#   per-sample projection) with the *integration variable* (which controls
#   what Harmony corrects for). This is both conceptually unsound and
#   unrealistic — in real multi-sample experiments, multiple samples typically
#   arise per batch.
#
#   To avoid this conflation, we generate N_PSEUDO_PER_BATCH (3) random
#   pseudo-samples within each batch, yielding 30 pseudo-samples total.
#   Cells within each batch are randomly partitioned with equal probability.
#   The pseudo-sample column is used only as `.sample.var`; Harmony still
#   corrects for `batch` via `.harmony.var`.
#
#   This design:
#     - Separates the sample-grouping axis from the batch-correction axis
#     - Tests tinydenseR's per-sample landmark selection at a realistic
#       granularity (~1,000–3,500 cells per pseudo-sample)
#     - Uses a fixed seed for reproducibility
#     - Does NOT affect Branch B (cell-level Harmony ignores sample structure)
###############################################################################

# Download: curl::curl_download("https://figshare.com/ndownloader/files/25717328", H5AD_PATH)

cat("\n── Preparing shared data ──\n")

# 1a. Read cell metadata from h5ad (handles both legacy and modern HDF5
#     categoricals via tinydenseR:::.h5ad_read_obs)
cell_meta     <- tinydenseR:::.h5ad_read_obs(H5AD_PATH)
gene_names    <- tinydenseR:::.h5ad_read_var_names(H5AD_PATH)
cell_barcodes <- tinydenseR:::.h5ad_read_obs_names(H5AD_PATH)

cat("  Cells:", nrow(cell_meta), " Genes:", length(gene_names), "\n")
cat("  Batches:", paste(sort(unique(as.character(cell_meta[[BATCH_VAR]]))),
                        collapse = ", "), "\n")
cat("  Cell types:", length(unique(cell_meta[[CELLTYPE_VAR]])), "\n")

# 1b. Open / cache BPCells on-disk matrix
if (dir.exists(BPCELLS_DIR)) {
  bp_mat <- tryCatch(BPCells::open_matrix_dir(BPCELLS_DIR), error = function(e) NULL)
}
if (!exists("bp_mat") || is.null(bp_mat)) {
  cat("  Converting h5ad matrix to BPCells on-disk format...\n")
  group <- tinydenseR:::.h5ad_resolve_counts_group(H5AD_PATH, NULL)
  h5_mat <- BPCells::open_matrix_anndata_hdf5(H5AD_PATH, group = group)
  bp_mat <- BPCells::write_matrix_dir(mat = h5_mat, dir = BPCELLS_DIR)
  bp_mat <- BPCells::open_matrix_dir(BPCELLS_DIR)
}
rownames(bp_mat) <- gene_names
colnames(bp_mat) <- cell_barcodes
cat("  BPCells matrix:", nrow(bp_mat), "×", ncol(bp_mat), "\n")

# 1c. Create pseudo-sample column
#     Within each batch, randomly assign cells to N_PSEUDO_PER_BATCH groups.
#     Names: "{batch}__ps{1..N}" (double-underscore avoids collision with batch
#     names that contain single underscores, e.g. "Sun_sample1_CS").
set.seed(SEED)
cell_meta[[PSEUDO_SAMPLE_VAR]] <- unsplit(
  lapply(split(seq_len(nrow(cell_meta)), cell_meta[[BATCH_VAR]]),
         function(idx) {
           batch_label <- cell_meta[[BATCH_VAR]][idx[1]]
           ps_ids <- sample(rep_len(seq_len(N_PSEUDO_PER_BATCH), length(idx)))
           paste0(batch_label, "__ps", ps_ids)
         }),
  cell_meta[[BATCH_VAR]]
)

n_pseudo <- length(unique(cell_meta[[PSEUDO_SAMPLE_VAR]]))
cat("  Pseudo-samples:", n_pseudo,
    "(", N_PSEUDO_PER_BATCH, "per batch ×", length(unique(cell_meta[[BATCH_VAR]])), "batches )\n")
cat("  Cells per pseudo-sample:\n")
print(table(cell_meta[[PSEUDO_SAMPLE_VAR]]))

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 2. BRANCH A: tinydenseR (landmark-level Harmony integration)
# ═══════════════════════════════════════════════════════════════════════════════
#
# We call RunTDR on the BPCells IterableMatrix with custom .cell.meta
# containing the pseudo_sample column. Key parameters:
#
#   .sample.var  = "pseudo_sample"  (30 pseudo-samples → landmark selection
#                                    and Symphony projection are per-pseudo-sample)
#   .harmony.var = "batch"          (Harmony corrects for the 10 real batches)
#
# This cleanly separates the sample-grouping axis from the batch-correction
# axis. The pipeline then:
#   a. Samples ~10% of cells as landmarks per pseudo-sample (leverage-score)
#   b. Computes PCA on landmarks (irlba, nPC=30, nHVG=5000)
#   c. Runs harmony::RunHarmony on landmark PCA embedding (vars_use="batch")
#   d. Builds symphony reference via symphony::buildReferenceFromHarmonyObj
#   e. Maps all cells via symphony::mapQuery (per pseudo-sample)
#   f. Computes UMAP for landmarks + transforms all cells
#
# After RunTDR, we extract the all-cell corrected PC embedding by replaying
# the same symphony::mapQuery that get.map uses internally, iterating over
# the per-pseudo-sample cell indices stored in tdr_obj@cells.
###############################################################################

cat("\n── Branch A: tinydenseR (landmark-level integration) ──\n")

set.seed(SEED)

tdr_obj <- tinydenseR::RunTDR(
  x               = bp_mat,
  .cell.meta      = cell_meta,
  .sample.var     = PSEUDO_SAMPLE_VAR,
  .assay.type     = "RNA",
  .harmony.var    = BATCH_VAR,
  .celltype.vec   = CELLTYPE_VAR,
  .prop.landmarks = PROP_LANDMARKS,
  .seed           = SEED,
  .verbose        = TRUE,
  .nHVG           = N_HVG,
  .nPC            = N_PC
)

# ── Extract Harmony settings used by tinydenseR ──────────────────────────────
# From landmarks.R ~line 859:
#   harmony::RunHarmony(
#     data_mat   = pca_embed,                  # landmark PC embedding
#     meta_data  = metadata[key, ],            # landmark-level metadata
#     vars_use   = harmony.var,                # = "batch"
#     return_obj = TRUE,
#     nclust     = min(round(nrow(pca_embed) / 30), 100) |> max(20),
#     verbose    = TRUE
#   )
# All other Harmony parameters are at their defaults:
#   theta = 2, lambda = 1, sigma = 0.1, max.iter.harmony = 10,
#   max.iter.cluster = 20, epsilon.harmony = 1e-4, etc.

n_landmarks <- nrow(tdr_obj@assay$expr)
tdr_nclust  <- min(round(n_landmarks / 30), 100) |> max(20)
cat("  Landmarks:", n_landmarks, "\n")
cat("  Harmony nclust (landmark-level):", tdr_nclust, "\n")
cat("  HVG count:", length(tdr_obj@landmark.embed$pca$HVG), "\n")
cat("  PCs:", ncol(tdr_obj@landmark.embed$pca$coord), "\n")

# ── Extract all-cell corrected embedding (TDR branch) ────────────────────────
# Replicate the per-(pseudo-)sample symphony::mapQuery projection that get.map
# uses (R/lm.graph.embed.R ~line 925). This gives the Harmony-corrected
# PC-space embedding for ALL cells, not just landmarks.

cat("  Projecting all cells through Symphony reference...\n")

harmony_ref <- tdr_obj@integration$harmony.obj
hvg_names   <- tdr_obj@landmark.embed$pca$HVG
pca_center  <- tdr_obj@landmark.embed$pca$center
pca_scale   <- tdr_obj@landmark.embed$pca$scale
raw_lm      <- tdr_obj@assay$raw   # landmarks raw counts (for size factor)

# Process each pseudo-sample exactly as get.map does.
# Track original cell indices directly (avoids fragile name-mangling).
embed_chunks   <- vector("list", length(tdr_obj@cells))
cell_idx_order <- integer(0)

for (i in seq_along(tdr_obj@cells)) {
  set.seed(SEED)
  
  col_idx <- tdr_obj@cells[[i]]
  
  # Extract raw counts for this pseudo-sample's cells
  mat_raw <- bp_mat[, col_idx, drop = FALSE]
  mat_raw <- methods::as(mat_raw, "dgCMatrix")
  mat_raw <- Matrix::t(mat_raw)
  
  n_cells <- 
    nrow(x = mat_raw)
  
  ## try rbindinf with landmarks
  #mat_raw <- rbind(mat_raw, raw_lm)

  # Size-factor normalization matching get.map (R/lm.graph.embed.R ~893-900)
  mat_norm <- (mat_raw / (Matrix::rowSums(mat_raw) /
                              mean(Matrix::rowSums(raw_lm))))
  
  mat_norm@x <- log2(mat_norm@x + 1)
  
  # create metadata for this pseudo-sample + landmarks
  #tmp_meta <-
  #  tdr_obj$metadata[tdr_obj$metadata$pseudo_sample == names(x = tdr_obj$cells)[[i]],] |>
  #  replicate(n = n_cells, simplify = FALSE) |>
  #  do.call(what = rbind) |>
  #  rbind(harmony_ref$meta_data)
  
  # Project through Symphony (same call as get.map: do_normalize=FALSE)
  Z <- symphony::mapQuery(
    exp_query       = Matrix::t(mat_norm[, hvg_names]),
    metadata_query  = matrix(data = 1, nrow = n_cells, ncol = 2),#tmp_meta
    ref_obj         = harmony_ref,
    vars            = NULL,#"batch",
    verbose         = FALSE,
    do_normalize    = FALSE,
    do_umap         = FALSE
  )$Z
  
  embed_chunks[[i]] <- Matrix::t(Z[,seq_len(length.out = n_cells)])
  cell_idx_order    <- c(cell_idx_order, col_idx)
}

tdr_allcell_embed <- do.call(rbind, embed_chunks)

cat("  TDR all-cell embedding:", nrow(tdr_allcell_embed), "×",
    ncol(tdr_allcell_embed), "\n")

# Map rows back to original cell_barcodes ordering via tracked indices
tdr_cell_ids <- cell_barcodes[cell_idx_order]

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 3. BRANCH B: Pure cell-level Harmony (matched settings)
# ═══════════════════════════════════════════════════════════════════════════════
#
# DESIGN DECISIONS for apples-to-apples comparison:
#
# MATCHED (identical to tinydenseR branch):
#   - Same h5ad file → same raw counts (via shared BPCells on-disk dir)
#   - Same HVG set (from tinydenseR's landmark-level selection)
#   - Same normalization: size-factor norm → log2(x + 1)
#   - Same PCA: irlba with same nPC=30, centering, scaling
#   - Same Harmony parameters: theta, lambda, sigma, max.iter, etc. (defaults)
#   - Same vars_use: "batch"
#   - Same seed
#
# NECESSARILY DIFFERENT:
#   - PCA is computed on ALL cells (not landmark subset)
#   - Harmony nclust is adapted to ALL cells: min(round(33506/30), 100) |> max(20)
#     = 100 (vs tinydenseR's ~tdr_nclust based on landmark count)
#   - Size-factor normalization: computed across all cells (no landmark
#     size-factor ratio adjustment, because there are no landmarks in this branch)
#
# RATIONALE for using tinydenseR's HVG set:
#   If we independently selected HVGs on all cells, the HVG sets would differ
#   (landmark-sampled population vs full population). This would conflate
#   HVG selection differences with integration method differences. By using
#   the SAME HVG set, any observed difference is attributable only to whether
#   integration operates on landmarks vs all cells.
#
# CAVEAT:
#   Using landmark-derived HVGs for the cell-level branch slightly favors
#   tinydenseR (its HVGs are tuned to its own landmarks). If this bias is
#   a concern, an alternative analysis could use all-cell HVGs for both —
#   but that would require rerunning tinydenseR with custom .force.in, which
#   could itself bias tinydenseR's landmark selection. No perfect solution
#   exists; we document the choice.
###############################################################################

cat("\n── Branch B: Cell-level Harmony (matched settings) ──\n")

# 3a. Materialize the count matrix for all cells on ALL genes first
#     (size-factor normalization must use total counts, not just HVG counts,
#     to match tinydenseR's normalization in landmarks.R)
cat("  Materializing full count matrix for normalization...\n")
mat_all_full <- methods::as(bp_mat, "dgCMatrix")
mat_all_full <- Matrix::t(mat_all_full)  # cells × genes
rownames(mat_all_full) <- cell_barcodes
colnames(mat_all_full) <- gene_names

# 3b. Normalize: same size-factor normalization as tinydenseR landmarks.R
#     tinydenseR computes size factors on ALL genes, then subsets to HVGs
cat("  Normalizing (size factors from all genes)...\n")
size_factors <- Matrix::rowSums(mat_all_full)
mat_norm_all <- mat_all_full / (size_factors / mean(size_factors))
mat_norm_all@x <- log2(mat_norm_all@x + 1)

# Now subset to HVGs (same order as tinydenseR)
cat("  Subsetting to", length(hvg_names), "HVGs...\n")
mat_norm_all <- mat_norm_all[, hvg_names]

# Free the full matrix
rm(mat_all_full); gc()

# 3c. PCA: same settings as tinydenseR landmarks.R pass 2
#     irlba, nv=30, center=colMeans, scale=colSds
cat("  Computing PCA (", N_PC, " PCs on all cells)...\n")

pca_center_all <- Matrix::colMeans(mat_norm_all)
pca_scale_all  <- sparseMatrixStats::colSds(mat_norm_all)

set.seed(SEED)
pca_all <- irlba::irlba(
  A      = mat_norm_all,
  nv     = N_PC,
  center = pca_center_all,
  scale  = pca_scale_all
)

pca_embed_all <- pca_all$u %*% diag(pca_all$d)
colnames(pca_embed_all) <- paste0("PC", seq_len(N_PC))
rownames(pca_embed_all) <- cell_barcodes

cat("  PCA embedding:", nrow(pca_embed_all), "×", ncol(pca_embed_all), "\n")

# 3d. Harmony: matched parameters
#     tinydenseR uses: nclust = min(round(n_landmarks/30), 100) |> max(20)
#     For all cells (33506): min(round(33506/30), 100) |> max(20) = 100
#     All other params are Harmony defaults (same as tinydenseR).
harmony_nclust_all <- min(round(nrow(pca_embed_all) / 30), 100) |> max(20)
cat("  Harmony nclust (cell-level):", harmony_nclust_all, "\n")
cat("  (tinydenseR used nclust =", tdr_nclust, "on", n_landmarks, "landmarks)\n")

set.seed(SEED)
harmony_all <- harmony::RunHarmony(
  data_mat      = pca_embed_all,
  meta_data     = cell_meta,
  vars_use      = BATCH_VAR,
  return_object = FALSE,    # returns corrected embedding matrix directly
  nclust        = harmony_nclust_all,
  verbose       = TRUE
)

colnames(harmony_all) <- paste0("PC", seq_len(N_PC))
rownames(harmony_all) <- cell_barcodes

cat("  Harmony cell-level embedding:", nrow(harmony_all), "×",
    ncol(harmony_all), "\n")

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 4. Align embeddings for comparison
# ═══════════════════════════════════════════════════════════════════════════════
#
# The TDR all-cell embedding was obtained by projecting each batch through
# the landmark-level Symphony reference. Cell ordering may differ from the
# original cell_meta ordering. We align here.
###############################################################################

cat("\n── Aligning embeddings ──\n")

# Ensure both embeddings cover the same cells
# tdr_cell_ids should match cell_barcodes (same cells, possibly different order)
# Some cells may have been dropped by tinydenseR due to min.cells.per.sample
tdr_coverage  <- sum(tdr_cell_ids %in% cell_barcodes)
cat("  TDR coverage:", tdr_coverage, "/", length(cell_barcodes), "cells\n")

if (tdr_coverage < length(cell_barcodes)) {
  warning(length(cell_barcodes) - tdr_coverage,
          " cells in the full dataset were not covered by tinydenseR. ",
          "Comparison will use the intersection.")
}

# Build aligned matrices (cells in cell_barcodes order, restricted to shared cells)
shared_cells <- intersect(cell_barcodes, tdr_cell_ids)
idx_tdr  <- match(shared_cells, tdr_cell_ids)
idx_harm <- match(shared_cells, cell_barcodes)

tdr_aligned   <- tdr_allcell_embed[idx_tdr, , drop = FALSE]
harm_aligned  <- harmony_all[idx_harm, , drop = FALSE]
meta_aligned  <- cell_meta[idx_harm, , drop = FALSE]

n_compared <- nrow(meta_aligned)
cat("  Cells in comparison:", n_compared, "\n")

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 5. Build shared kNN graphs for metric computation
# ═══════════════════════════════════════════════════════════════════════════════
###############################################################################

cat("\n── Building kNN graphs ──\n")

K_NN <- 30L  # k for nearest-neighbor graph

# TDR branch kNN (using RANN for exact kNN; k+1 because nn2 includes self)
knn_tdr_raw <- RANN::nn2(
  data  = as.matrix(tdr_aligned),
  query = as.matrix(tdr_aligned),
  k     = K_NN + 1L
)
# Remove self-match (first column)
knn_tdr <- list(idx = knn_tdr_raw$nn.idx[, -1, drop = FALSE],
                dist = knn_tdr_raw$nn.dists[, -1, drop = FALSE])

# Harmony branch kNN
knn_harm_raw <- RANN::nn2(
  data  = as.matrix(harm_aligned),
  query = as.matrix(harm_aligned),
  k     = K_NN + 1L
)
knn_harm <- list(idx = knn_harm_raw$nn.idx[, -1, drop = FALSE],
                 dist = knn_harm_raw$nn.dists[, -1, drop = FALSE])

cat("  kNN built for both branches (k =", K_NN, ")\n")

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 6. Integration quality metrics
# ═══════════════════════════════════════════════════════════════════════════════
#
# METRIC RATIONALE:
#
# We compute MCC + 5 additional metrics addressing both:
#   (A) Biological conservation — does integration preserve cell-type structure?
#   (B) Batch mixing — does integration remove unwanted batch effects?
#
# ┌───────────────────────────┬───────────┬─────────────────────────────────────┐
# │ Metric                    │ Category  │ What it captures                    │
# ├───────────────────────────┼───────────┼─────────────────────────────────────┤
# │ 1. MCC (annotation)       │ Bio cons. │ Label agreement after kNN transfer  │
# │ 2. Silhouette (celltype)  │ Bio cons. │ Compactness of cell types           │
# │ 3. Silhouette (batch)     │ Batch mix │ Mixing of batches (lower = better)  │
# │ 4. kBET acceptance rate   │ Batch mix │ Local batch composition vs expected │
# │ 5. Graph connectivity     │ Bio cons. │ Connectedness of cell-type clusters │
# │ 6. Batch entropy of kNN   │ Batch mix │ Diversity of batches in local kNN   │
# └───────────────────────────┴───────────┴─────────────────────────────────────┘
#
# Metric details:
#
# 1. MCC (Matthews Correlation Coefficient) on kNN-transferred labels:
#    For each cell, assign a label by majority vote from its k nearest
#    neighbors (leaving one out). Compare to true label. MCC is robust to
#    class imbalance and ranges [-1, 1] (1 = perfect). Higher = better
#    biological conservation.
#
# 2. Average Silhouette Width (ASW) by cell type:
#    Measures how well cell types form distinct clusters in the embedding.
#    Range [-1, 1]; higher = tighter, better-separated cell types.
#    Captures whether integration preserves major biological axes.
#
# 3. Average Silhouette Width (ASW) by batch (1 - |ASW_batch|):
#    Measures batch mixing. Perfect mixing → ASW_batch ≈ 0 → score ≈ 1.
#    Strong batch effect → |ASW_batch| high → score low.
#    Rescaled to [0, 1] where higher = better batch mixing.
#
# 4. kBET acceptance rate:
#    Tests whether local neighborhoods have the same batch composition
#    as the global dataset (chi-squared test). Reports the fraction of
#    cells that "pass" (accept H0 that local composition = global).
#    Higher = better batch mixing. Specifically sensitive to local
#    mixing quality, complementing the global ASW_batch metric.
#
# 5. Graph connectivity score:
#    For each cell type, compute the fraction of cells in the largest
#    connected component of the cell-type-specific kNN subgraph.
#    Average across cell types. Higher = better preservation of cell-type
#    topology. Captures fragmentation that silhouette misses.
#
# 6. Batch entropy of mixing (Shannon):
#    For each cell, compute the Shannon entropy of batch labels among
#    its k nearest neighbors, normalized by max possible entropy.
#    Average across cells. Range [0, 1]; higher = more uniform batch
#    mixing in local neighborhoods. Complements kBET by being
#    non-parametric and not dependent on chi-squared assumptions.
###############################################################################

cat("\n── Computing integration metrics ──\n")

# ── Helper functions ──────────────────────────────────────────────────────────

#' kNN label transfer and MCC computation
compute_knn_mcc <- function(knn_idx, true_labels) {
  # For each cell, predict label by majority vote among its k nearest neighbors
  neighbor_labels <- matrix(true_labels[knn_idx], nrow = nrow(knn_idx), ncol = ncol(knn_idx))
  pred <- apply(neighbor_labels, 1, function(row) {
    tbl <- table(row)
    names(tbl)[which.max(tbl)]
  })
  # Multiclass MCC (sklearn-style via confusion matrix)
  .multiclass_mcc(true_labels, pred)
}

#' Multiclass MCC via the confusion matrix formulation
.multiclass_mcc <- function(y_true, y_pred) {
  cm <- table(True = y_true, Pred = y_pred)
  n  <- sum(cm)
  # Sum of correct predictions
  c_sum <- sum(diag(cm))
  # Row and column marginals
  s  <- sum(cm^2)
  p  <- colSums(cm)
  t_ <- rowSums(cm)
  
  numerator   <- c_sum * n - sum(p * t_)
  denominator <- sqrt((n^2 - sum(p^2)) * (n^2 - sum(t_^2)))
  
  if (denominator == 0) return(0)
  numerator / denominator
}

#' Average silhouette width (using Euclidean distance on embedding)
#' Subsampled for computational efficiency
compute_asw <- function(embed, labels, n_subsample = 5000L) {
  set.seed(SEED)
  if (nrow(embed) > n_subsample) {
    idx <- sample(nrow(embed), n_subsample)
    embed  <- embed[idx, ]
    labels <- labels[idx]
  }
  d <- dist(embed)
  sil <- cluster::silhouette(as.integer(factor(labels)), d)
  mean(sil[, "sil_width"])
}

#' kBET acceptance rate
#' Tests whether local neighborhoods have the expected batch composition
compute_kbet <- function(knn_idx, batch_labels, n_subsample = 1000L, k0 = NULL) {
  if (is.null(k0)) k0 <- min(25L, ncol(knn_idx))
  
  set.seed(SEED)
  n <- length(batch_labels)
  idx_sub <- sample(n, min(n_subsample, n))
  
  # Global batch frequencies
  global_freq <- table(batch_labels) / n
  batches <- names(global_freq)
  
  n_accept <- 0L
  for (i in idx_sub) {
    neighbors <- knn_idx[i, seq_len(k0)]
    local_counts <- table(factor(batch_labels[neighbors], levels = batches))
    # Chi-squared test against expected frequencies
    expected <- global_freq * k0
    # Guard against zero expected
    valid <- expected > 0
    if (sum(valid) < 2) {
      n_accept <- n_accept + 1L
      next
    }
    chi2 <- sum((local_counts[valid] - expected[valid])^2 / expected[valid])
    p_val <- stats::pchisq(chi2, df = sum(valid) - 1, lower.tail = FALSE)
    if (p_val > 0.05) n_accept <- n_accept + 1L
  }
  
  n_accept / length(idx_sub)
}

#' Graph connectivity: fraction of cells in largest connected component
#' per cell type, averaged across types
compute_graph_connectivity <- function(knn_idx, celltype_labels) {
  celltypes <- unique(celltype_labels)
  scores <- vapply(celltypes, function(ct) {
    ct_mask <- which(celltype_labels == ct)
    if (length(ct_mask) < 2) return(1.0)
    
    # Build adjacency among ct_mask cells only
    # Map global → local indices
    global_to_local <- rep(NA_integer_, length(celltype_labels))
    global_to_local[ct_mask] <- seq_along(ct_mask)
    
    n_local <- length(ct_mask)
    # Build edge list
    edges_from <- integer(0)
    edges_to   <- integer(0)
    for (j in seq_along(ct_mask)) {
      gi <- ct_mask[j]
      neighbors <- knn_idx[gi, ]
      local_neighbors <- global_to_local[neighbors]
      local_neighbors <- local_neighbors[!is.na(local_neighbors)]
      if (length(local_neighbors) > 0) {
        edges_from <- c(edges_from, rep(j, length(local_neighbors)))
        edges_to   <- c(edges_to, local_neighbors)
      }
    }
    
    if (length(edges_from) == 0) return(0.0)
    
    g <- igraph::graph_from_edgelist(
      cbind(edges_from, edges_to),
      directed = FALSE
    )
    # Ensure all ct cells are represented as vertices (some may be isolated)
    g <- igraph::add_vertices(g, max(0, n_local - igraph::vcount(g)))
    g <- igraph::simplify(g)
    
    comps <- igraph::components(g)
    max(comps$csize) / n_local
  }, numeric(1))
  
  mean(scores)
}

#' Batch entropy of mixing (normalized Shannon entropy of batch in kNN)
compute_batch_entropy <- function(knn_idx, batch_labels) {
  batches <- unique(batch_labels)
  n_batches <- length(batches)
  max_entropy <- log(n_batches)
  
  if (max_entropy == 0) return(1.0)
  
  n <- nrow(knn_idx)
  entropies <- vapply(seq_len(n), function(i) {
    neighbor_batches <- batch_labels[knn_idx[i, ]]
    freq <- table(factor(neighbor_batches, levels = batches)) / length(neighbor_batches)
    freq <- freq[freq > 0]
    -sum(freq * log(freq)) / max_entropy
  }, numeric(1))
  
  mean(entropies)
}

# ── Compute all metrics ──────────────────────────────────────────────────────

true_ct    <- meta_aligned[[CELLTYPE_VAR]]
true_batch <- meta_aligned[[BATCH_VAR]]

cat("  [1/6] MCC (kNN label transfer)...\n")
mcc_tdr  <- compute_knn_mcc(knn_tdr$idx, true_ct)
mcc_harm <- compute_knn_mcc(knn_harm$idx, true_ct)

cat("  [2/6] Silhouette width (cell type)...\n")
asw_ct_tdr  <- compute_asw(as.matrix(tdr_aligned), true_ct)
asw_ct_harm <- compute_asw(as.matrix(harm_aligned), true_ct)

cat("  [3/6] Silhouette width (batch, inverted)...\n")
asw_batch_tdr  <- 1 - abs(compute_asw(as.matrix(tdr_aligned), true_batch))
asw_batch_harm <- 1 - abs(compute_asw(as.matrix(harm_aligned), true_batch))

cat("  [4/6] kBET acceptance rate...\n")
kbet_tdr  <- compute_kbet(knn_tdr$idx, true_batch)
kbet_harm <- compute_kbet(knn_harm$idx, true_batch)

cat("  [5/6] Graph connectivity...\n")
gc_tdr  <- compute_graph_connectivity(knn_tdr$idx, true_ct)
gc_harm <- compute_graph_connectivity(knn_harm$idx, true_ct)

cat("  [6/6] Batch entropy of mixing...\n")
be_tdr  <- compute_batch_entropy(knn_tdr$idx, true_batch)
be_harm <- compute_batch_entropy(knn_harm$idx, true_batch)

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 7. Results summary
# ═══════════════════════════════════════════════════════════════════════════════
###############################################################################

cat("\n── Results ──\n\n")

results_df <- data.frame(
  Metric = c(
    "MCC (cell type, kNN transfer)",
    "ASW cell type (bio conservation)",
    "ASW batch (1 - |ASW|, batch mixing)",
    "kBET acceptance rate (batch mixing)",
    "Graph connectivity (bio conservation)",
    "Batch entropy of mixing"
  ),
  Category = c("Bio", "Bio", "Batch", "Batch", "Bio", "Batch"),
  Higher_is_better = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  tinydenseR = round(c(mcc_tdr, asw_ct_tdr, asw_batch_tdr, kbet_tdr, gc_tdr, be_tdr), 4),
  Harmony_cellevel = round(c(mcc_harm, asw_ct_harm, asw_batch_harm, kbet_harm, gc_harm, be_harm), 4),
  stringsAsFactors = FALSE
)

results_df$Delta <- results_df$tinydenseR - results_df$Harmony_cellevel

print(results_df, row.names = FALSE)

# Save results
write.csv(results_df,
          file = file.path(RESULTS_DIR, "integration_metrics.csv"),
          row.names = FALSE)

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 8. Visualization
# ═══════════════════════════════════════════════════════════════════════════════
###############################################################################

cat("\n── Generating visualizations ──\n")

# 8a. UMAP for visual comparison
cat("  Computing UMAPs...\n")
set.seed(SEED)
umap_tdr <- uwot::umap(
  X         = as.matrix(tdr_aligned),
  n_threads = N_THREADS,
  seed      = SEED,
  verbose   = FALSE
)
rownames(umap_tdr) <- rownames(tdr_aligned)

set.seed(SEED)
umap_harm <- uwot::umap(
  X         = as.matrix(harm_aligned),
  n_threads = N_THREADS,
  seed      = SEED,
  verbose   = FALSE
)
rownames(umap_harm) <- rownames(harm_aligned)

# Build plot data
plot_df <- rbind(
  data.frame(
    UMAP1 = umap_tdr[, 1], UMAP2 = umap_tdr[, 2],
    cell_type = true_ct, batch = true_batch,
    method = "tinydenseR (landmarks)"
  ),
  data.frame(
    UMAP1 = umap_harm[, 1], UMAP2 = umap_harm[, 2],
    cell_type = true_ct, batch = true_batch,
    method = "Harmony (all cells)"
  )
)

# 8b. UMAP colored by cell type
p_ct <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = I(x = 0.1)) +
  facet_wrap(~ method) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(ncol = 4, override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Cell type preservation",
       subtitle = "Luecken et al. Immune_ALL_human (33,506 cells, 10 batches)") +
  ggh4x::force_panelsizes(cols = unit(3, "in"),
                          rows = unit(3, "in"))

# 8c. UMAP colored by batch
p_batch <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = I(x = 0.1)) +
  facet_wrap(~ method) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(ncol = 5, override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Batch mixing",
       subtitle = "Luecken et al. Immune_ALL_human (33,506 cells, 10 batches)") +
  ggh4x::force_panelsizes(cols = unit(3, "in"),
                          rows = unit(3, "in"))

# 8d. Metric comparison bar chart
metric_long <- results_df |>
  tidyr::pivot_longer(
    cols = c(tinydenseR, Harmony_cellevel),
    names_to = "Method",
    values_to = "Score"
  ) |>
  dplyr::mutate(
    Method = dplyr::recode(Method,
                           tinydenseR = "tinydenseR (landmarks)",
                           Harmony_cellevel = "Harmony (all cells)")
  )

p_metrics <- ggplot(metric_long, aes(x = Metric, y = Score, fill = Method)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("tinydenseR (landmarks)" = "#2166AC",
                                "Harmony (all cells)" = "#B2182B")) +
  labs(title = "Integration quality metrics",
       subtitle = "Higher is better for all metrics",
       y = "Score", x = NULL)  +
  ggh4x::force_panelsizes(cols = unit(3, "in"),
                          rows = unit(1.5, "in"))

# Combine
p_combined <- (p_ct / p_batch / p_metrics) +
  plot_annotation(
    title = "Landmark-level vs cell-level integration",
    subtitle = paste("tinydenseR RunTDR (", n_landmarks, "landmarks, Harmony on landmarks)\n",
                     "vs Harmony on all", n_compared, "cells"),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(filename = file.path(RESULTS_DIR, "integration_comparison.pdf"),
       plot = p_combined,
       width = 10, height = 15, units = "in", dpi = 300)

ggsave(filename = file.path(RESULTS_DIR, "integration_comparison.pdf"),
       plot = p_combined,
       width = 10, height = 15, units = "in", dpi = 300)

(p <-
  p_ct + 
  labs(
    color = "final_annotation") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank()) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(size = I(x = 5)))) +
  facet_grid(
    cols = vars(method),
    labeller = as_labeller(
      x =
         c(
          "tinydenseR (landmarks)" = "landmark-level",
          "Harmony (all cells)" = "cell-level")))  +
  force_panelsizes(
    cols = unit(2, "in"),
    rows = unit(2, "in")) +
  scale_color_manual(
    values = colorRampPalette(
      colors = unname(
        obj = tinydenseR::Color.Palette[1,1:5]))(
          length(x = unique(x = p_ct$data$cell_type))))); ggsave(
    plot = p,
    filename = file.path(RESULTS_DIR, "p_ct.png"),
    width = 9,
    height = 3,
    units = "in",
    dpi = 300); rm(p)

(p <-
  p_batch + 
  labs(
    color = "batch") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank()) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(size = I(x = 5)))) +
  facet_grid(
    cols = vars(method),
    labeller = as_labeller(
      x =
         c(
          "tinydenseR (landmarks)" = "landmark-level",
          "Harmony (all cells)" = "cell-level")))  +
  force_panelsizes(
    cols = unit(2, "in"),
    rows = unit(2, "in")) +
  scale_color_manual(
    values = colorRampPalette(
      colors = unname(
        obj = tinydenseR::Color.Palette[1,1:5]))(
          length(x = unique(x = p_batch$data$cell_type))))); ggsave(
    plot = p,
    filename = file.path(RESULTS_DIR, "p_batch.png"),
    width = 7.5,
    height = 3,
    units = "in",
    dpi = 300); rm(p)

(p <-
  p_metrics +
  labs(
    title = "Integration quality metrics") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_blank()) +
  scale_fill_manual(
    values = unname(obj = tinydenseR::Color.Palette[1,1:2]),
    labels = c(
          "tinydenseR (landmarks)" = "landmark-level",
          "Harmony (all cells)" = "cell-level"))); ggsave(
    plot = p,
    filename = file.path(RESULTS_DIR, "p_metrics.png"),
    width = 7,
    height = 2.5,
    units = "in",
    dpi = 300); rm(p)

cat("  Saved to:", RESULTS_DIR, "\n")

###############################################################################
# ═══════════════════════════════════════════════════════════════════════════════
# 9. Session info & comparison design documentation
# ═══════════════════════════════════════════════════════════════════════════════
###############################################################################

cat("\n── Comparison design summary ──\n\n")
cat("MATCHED between branches:\n")
cat("  - Raw count matrix: identical (BPCells on-disk)\n")
cat("  - HVG set:", length(hvg_names), "genes (from tinydenseR landmark selection)\n")
cat("  - Normalization: size-factor + log2(x+1)\n")
cat("  - PCA method: irlba, nPC =", N_PC, "\n")
cat("  - Harmony vars_use:", BATCH_VAR, "\n")
cat("  - Harmony defaults: theta=2, lambda=1, sigma=0.1 (both branches)\n")
cat("  - Random seed:", SEED, "\n")

cat("\nNECESSARILY DIFFERENT:\n")
cat("  - PCA computed on:", n_landmarks, "landmarks (TDR) vs", n_compared, "cells (Harmony)\n")
cat("  - Harmony inputs:", n_landmarks, "landmark PCs (TDR) vs", n_compared, "cell PCs (Harmony)\n")
cat("  - Harmony nclust:", tdr_nclust, "(TDR) vs", harmony_nclust_all, "(Harmony)\n")
cat("  - TDR all-cell embedding obtained via Symphony projection (approximate)\n")
cat("  - Size-factor normalization differs slightly (landmark vs global mean)\n")

cat("\nSAMPLE STRUCTURE:\n")
cat("  - tinydenseR .sample.var:", PSEUDO_SAMPLE_VAR, "(", n_pseudo, "pseudo-samples,",
    N_PSEUDO_PER_BATCH, "per batch )\n")
cat("  - tinydenseR .harmony.var:", BATCH_VAR, "(", length(unique(cell_meta[[BATCH_VAR]])), "real batches )\n")
cat("  - Pseudo-samples separate the grouping axis (landmark selection, Symphony\n")
cat("    projection) from the correction axis (Harmony batch integration).\n")
cat("  - Cell-level Harmony (Branch B) does not use sample structure.\n")

cat("\nLIMITATION:\n")
cat("  The tinydenseR all-cell embedding is obtained by projecting through the\n")
cat("  landmark-level Symphony reference. This is the SAME operation get.map()\n")
cat("  uses internally. It is an approximation: the correction model was fit on\n")
cat("  landmarks only and applied to all cells via linear projection. This is\n")
cat("  inherent to the landmark-based design and IS the feature being evaluated.\n")

cat("\n── Session info ──\n")
sessionInfo()

cat("\n── Done ──\n")
