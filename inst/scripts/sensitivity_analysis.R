###############################################################################
# Hyperparameter Sensitivity Analysis for tinydenseR
#
# Purpose:
#   Systematically evaluate the stability of tinydenseR results across key
#   hyperparameter choices on the COVID-19 COMBAT dataset. One parameter is
#   varied at a time while all others are held at their defaults.
#
# Design:
#   - Parameter grid: .tot.landmarks, .prop.landmarks, .k, .nHVG, .nPC
#   - Response variables: density contrast (Spearman rho vs default),
#     pePC1 variance explained and score correlation, plsD1 score and
#     loading correlations, wall-clock time.
#   - Fixed model: early CRIT vs early MILD (one contrast per run)
#   - All runs use seed = 123 for reproducibility
#
# Usage:
#   Rscript inst/scripts/sensitivity_analysis.R \
#     --h5ad ~/Downloads/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad \
#     --outdir results/sensitivity
#
# Output:
#   An RDS file containing a list of tidy data.frames suitable for plotting.
#
# Requirements:
#   tinydenseR, anndataR, BPCells, dplyr, limma
#
# Copyright 2025 Novartis Biomedical Research Inc.
# SPDX-License-Identifier: MIT
###############################################################################

# ===========================================================================
# 0. Parse command-line arguments
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)

# Defaults
h5ad_path <- "~/Downloads/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad"
out_dir   <- "results/sensitivity"

if ("--h5ad" %in% args) {
  h5ad_path <- args[which(args == "--h5ad") + 1L]
}
if ("--outdir" %in% args) {

  out_dir <- args[which(args == "--outdir") + 1L]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load packages
# ===========================================================================

library(tinydenseR)
library(dplyr)
library(limma)
library(anndataR)
library(BPCells)

# ===========================================================================
# 2. Define parameter grid (one-at-a-time design)
# ===========================================================================

# Default values (matching vignette / package defaults)
defaults <- list(
  .tot.landmarks  = 5000,
  .prop.landmarks = 0.10,

  .k              = 20,
  .nHVG           = 5000,
  .nPC            = 30
)

# Parameter sweep values
param_grid <- list(
  .tot.landmarks  = c(1000, 2500, 5000, 7500, 10000),
  .prop.landmarks = c(0.01, 0.05, 0.10, 0.15, 0.20),
  .k              = c(10, 15, 20, 30, 50),
  .nHVG           = c(1000, 2000, 5000, 8000),
  .nPC            = c(10, 15, 20, 30, 50)
)

# Build run table: one row per (parameter, value) pair
run_table <- do.call(rbind, lapply(names(param_grid), function(param) {
  data.frame(
    parameter = param,
    value     = param_grid[[param]],
    is_default = param_grid[[param]] == defaults[[param]],
    stringsAsFactors = FALSE
  )
}))

message(sprintf("Total runs: %d (including %d default-value runs)",
                nrow(run_table),
                sum(run_table$is_default)))

# ===========================================================================
# 3. Load and prepare data (same preprocessing as vignette)
# ===========================================================================

message("Loading COMBAT data...")

h5ad <- anndataR::read_h5ad(path = h5ad_path, as = "HDF5AnnData")
feature_metadata <- h5ad$var
metadata <- h5ad$obs

# Basic QC: minimum 1000 UMI and 500 genes per cell
keep_QC_cells <-
  (metadata$QC_total_UMI >= 1000) &
  (metadata$QC_ngenes >= 500)

# Keep only HV, MILD, SEV, and CRIT patients with >= 200 cells
keep.pt.group <-
  metadata$scRNASeq_sample_ID[
    metadata$Source %in% c("HV", "COVID_MILD", "COVID_SEV", "COVID_CRIT")
  ] |> 
  table() |>
  (\(x) names(x[x >= 200]))()

smpl.to.keep <-
  keep.QC.cells & (metadata$scRNASeq_sample_ID %in% keep.pt.group)

metadata <-
  metadata[smpl.to.keep,]

# Encode covariates
metadata <- metadata |>
  dplyr::mutate(
    Source = factor(
      x = Source,
      levels = c("HV", "COVID_MILD", "COVID_SEV", "COVID_CRIT")),
    TimeSinceOnset = ifelse(
      test = Source != "HV",
      yes = TimeSinceOnset,
      no = -1)
  ) |>
  droplevels()

# On-disk count matrix via BPCells
keep_unique_genes <- !duplicated(feature_metadata$feature_name)
feature_metadata <- feature_metadata[keep_unique_genes, ] |> droplevels()

ondisk_path <-
  dirname(path = h5ad_path) |>
  file.path("COMBAT_COVID_ondisk")

ondisk_data <-
  BPCells::open_matrix_anndata_hdf5(path = h5ad_path, group = "/X")

if (!file.exists(ondisk_path)) {
  BPCells::write_matrix_dir(mat = ondisk_data, dir = ondisk_path)
}

ondisk_mat <- BPCells::open_matrix_dir(dir = ondisk_path)

ondisk.mat <-
  ondisk.mat[keep.unique.genes, smpl.to.keep]

rownames(ondisk_mat) <- feature_metadata$feature_name

message("Data loaded: ", ncol(ondisk_mat), " cells, ", nrow(ondisk_mat), " genes")

# ===========================================================================
# 4. Helper: run one configuration and extract metrics
# ===========================================================================

#' Run a single tinydenseR configuration and return metrics
#'
#' @param params Named list of parameters to override (e.g., list(.k = 30))
#' @param ondisk_mat On-disk count matrix
#' @param metadata Cell-level metadata data.frame
#' @return Named list with metrics or NULL on failure
run_one_config <- function(params, ondisk_mat, metadata) {

  t_start <- proc.time()["elapsed"]

  # Merge with defaults
  run_params <- defaults
  run_params[names(params)] <- params

  # --- Step 1: RunTDR ---
  lm_cells <- tryCatch(
    tinydenseR::RunTDR(
      x = ondisk_mat,
      .sample.var      = "scRNASeq_sample_ID",
      .cell.meta       = metadata,
      .assay.type      = "RNA",
      .celltype.vec    = "cell_type",
      .label.confidence = 0.5,
      .verbose         = FALSE,
      .seed            = 123,
      .tot.landmarks   = run_params$.tot.landmarks,
      .prop.landmarks  = run_params$.prop.landmarks,
      .k               = run_params$.k,
      .nHVG            = run_params$.nHVG,
      .nPC             = run_params$.nPC
    ),
    error = function(e) {
      message("  RunTDR failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(lm_cells)) return(NULL)

  # --- Step 2: Prepare covariates and design (same as vignette) ---
  lm_cells@metadata <- dplyr::mutate(
    .data = lm_cells@metadata,
    TSO.binary = dplyr::case_when(
      Source == "HV" ~ "uninfected",
      Source != "HV" ~ ifelse(
        TimeSinceOnset > median(TimeSinceOnset[Source != "HV"]),
        "late", "early")
    ) |> factor(levels = c("uninfected", "early", "late"))
  )

  lm_cells@metadata$Age <-
    make.names(lm_cells@metadata$Age) |>
    gsub(pattern = "^X", replacement = "Age.", fixed = FALSE)

  lm_cells@metadata$Source_TSO.binary <-
    paste(lm_cells@metadata$Source,
          lm_cells@metadata$TSO.binary, sep = "_") |>
    gsub(pattern = "^HV_|^COVID_", replacement = "", fixed = FALSE) |>
    factor(levels = c("uninfected", "MILD_early", "MILD_late",
                      "SEV_early", "SEV_late", "CRIT_early", "CRIT_late"))

  .design <- tryCatch(
    model.matrix(~ 0 + Source_TSO.binary + sex + Age,
                 data = lm_cells@metadata) |>
      (\(x) {
        colnames(x) <- gsub("^Source_TSO.binary", "", colnames(x))
        x
      })(),
    error = function(e) {
      message("  Design matrix construction failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(.design)) return(NULL)

  .contrast.matrix <- tryCatch(
    limma::makeContrasts(
      CRIT_earlyVSMILD_early = CRIT_early - MILD_early,
      levels = .design
    ),
    error = function(e) {
      message("  Contrast matrix failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(.contrast.matrix)) return(NULL)

  # --- Step 3: get.lm ---
  lm_cells <- tryCatch(
    tinydenseR::get.lm(
      x = lm_cells,
      .design = .design,
      .contrasts = .contrast.matrix,
      .model.name = "sensitivity",
      .verbose = FALSE
    ),
    error = function(e) {
      message("  get.lm failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(lm_cells)) return(NULL)

  # --- Step 4: get.embedding (pePC) ---
  lm_cells <- tryCatch(
    tinydenseR::get.embedding(
      x = lm_cells,
      .full.model = "sensitivity",
      .contrast.of.interest = "CRIT_earlyVSMILD_early",
      .verbose = FALSE
    ),
    error = function(e) {
      message("  get.embedding failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(lm_cells)) return(NULL)

  # --- Step 5: get.plsD ---
  lm_cells <- tryCatch(
    tinydenseR::get.plsD(
      x = lm_cells,
      .coef.col = "CRIT_earlyVSMILD_early",
      .model.name = "sensitivity",
      .verbose = FALSE,
      .min.prop = 0.001
    ),
    error = function(e) {
      message("  get.plsD failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(lm_cells)) return(NULL)

  t_elapsed <- proc.time()["elapsed"] - t_start

  # --- Extract metrics ---
  # Density contrast vector (log fold-change per landmark)
  density_contrast <-
    lm_cells@results$lm$sensitivity$fit$coefficients[, "CRIT_earlyVSMILD_early"]

  # pePC1 scores and variance explained
  pepc_slot <- lm_cells@sample.embed$pepc[["CRIT_earlyVSMILD_early"]]
  pePC1_scores <- pepc_slot$coord[, "pePC1", drop = TRUE]  # named vector (samples)
  pePC1_var_explained <- pepc_slot$perc.tot.var.exp["pePC1"]

  # plsD1 scores and loadings
  pls_slot <- lm_cells@results$pls[["CRIT_earlyVSMILD_early"]]
  plsD1_scores <- pls_slot$coord[, "plsD1"]
  plsD1_loadings <- pls_slot$loadings[, "plsD1"]

  list(
    density_contrast    = density_contrast,
    pePC1_scores        = pePC1_scores,
    pePC1_var_explained = pePC1_var_explained,
    plsD1_scores        = plsD1_scores,
    plsD1_loadings      = plsD1_loadings,
    elapsed_sec         = unname(t_elapsed)
  )
}

# ===========================================================================
# 5. Run the default configuration first (reference baseline)
# ===========================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Running DEFAULT configuration...")
message(paste(rep("=", 70), collapse = ""))

set.seed(123) # no-op safety measure
ref_result <- run_one_config(params = list(), ondisk_mat = ondisk_mat, metadata = metadata)

if (is.null(ref_result)) {

  stop("Default configuration failed. Cannot proceed with sensitivity analysis.")
}

message(sprintf("Default run completed in %.1f seconds", ref_result$elapsed_sec))

# ===========================================================================
# 6. Main sweep loop
# ===========================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Starting parameter sweep...")
message(paste(rep("=", 70), collapse = ""))

results_list <- vector("list", nrow(run_table))

for (i in seq_len(nrow(run_table))) {

  param_name  <- run_table$parameter[i]
  param_value <- run_table$value[i]

  message(sprintf("\n[%d/%d] %s = %s %s",
                  i, nrow(run_table), param_name, param_value,
                  if (run_table$is_default[i]) "(default)" else ""))

  # Skip actual computation for default value (reuse ref_result)
  if (run_table$is_default[i]) {
    run_result <- ref_result
  } else {
    params_override <- stats::setNames(list(param_value), param_name)
    set.seed(123)
    run_result <- run_one_config(
      params     = params_override,
      ondisk_mat = ondisk_mat,
      metadata   = metadata
    )
  }

  if (is.null(run_result)) {
    results_list[[i]] <- data.frame(
      parameter             = param_name,
      value                 = param_value,
      is_default            = run_table$is_default[i],
      density_rho           = NA_real_,
      pePC1_var_explained   = NA_real_,
      pePC1_score_rho       = NA_real_,
      plsD1_score_rho       = NA_real_,
      plsD1_loading_rho     = NA_real_,
      elapsed_sec           = NA_real_,
      status                = "failed",
      stringsAsFactors      = FALSE
    )
    next
  }

  # --- Compute Spearman correlations vs. reference ---


  # Density contrast: landmark-level correlation is only meaningful when the

  # landmark set is identical (same size AND same sampling seed). When
  # .tot.landmarks differs from the default, landmarks are resampled at a

  # different resolution so `intersect()` yields a biased subset — report NA
  # and rely on the sample-level pePC metric instead.
  landmarks_comparable <-
    !(param_name == ".tot.landmarks" &&
        param_value != defaults$.tot.landmarks)

  if (landmarks_comparable) {
    shared_lm <- intersect(names(run_result$density_contrast),
                           names(ref_result$density_contrast))
    if (length(shared_lm) > 10 &&
        length(shared_lm) == length(ref_result$density_contrast)) {
      density_rho <- cor(
        run_result$density_contrast[shared_lm],
        ref_result$density_contrast[shared_lm],
        method = "spearman"
      )
    } else {
      density_rho <- NA_real_
    }
  } else {
    density_rho <- NA_real_
  }

  # pePC1 scores: align by sample names
  shared_samples <- intersect(names(run_result$pePC1_scores),
                              names(ref_result$pePC1_scores))
  if (length(shared_samples) > 3) {
    pePC1_score_rho <- cor(
      run_result$pePC1_scores[shared_samples],
      ref_result$pePC1_scores[shared_samples],
      method = "spearman"
    )
  } else {
    pePC1_score_rho <- NA_real_
  }

  # plsD1 scores: per-landmark quantity — same comparability constraint
  if (landmarks_comparable) {
    shared_lm_pls <- intersect(names(run_result$plsD1_scores),
                               names(ref_result$plsD1_scores))
    if (length(shared_lm_pls) > 10 &&
        length(shared_lm_pls) == length(ref_result$plsD1_scores)) {
      plsD1_score_rho <- cor(
        run_result$plsD1_scores[shared_lm_pls],
        ref_result$plsD1_scores[shared_lm_pls],
        method = "spearman"
      )
    } else {
      plsD1_score_rho <- NA_real_
    }
  } else {
    plsD1_score_rho <- NA_real_
  }

  # plsD1 loadings: align by gene names
  shared_genes <- intersect(names(run_result$plsD1_loadings),
                            names(ref_result$plsD1_loadings))
  if (length(shared_genes) > 10) {
    plsD1_loading_rho <- cor(
      run_result$plsD1_loadings[shared_genes],
      ref_result$plsD1_loadings[shared_genes],
      method = "spearman"
    )
  } else {
    plsD1_loading_rho <- NA_real_
  }

  results_list[[i]] <- data.frame(
    parameter             = param_name,
    value                 = param_value,
    is_default            = run_table$is_default[i],
    density_rho           = density_rho,
    pePC1_var_explained   = run_result$pePC1_var_explained,
    pePC1_score_rho       = pePC1_score_rho,
    plsD1_score_rho       = plsD1_score_rho,
    plsD1_loading_rho     = plsD1_loading_rho,
    elapsed_sec           = run_result$elapsed_sec,
    status                = "success",
    stringsAsFactors      = FALSE
  )

  message(sprintf("  density_rho=%.3f | pePC1_var=%.1f%% | pePC1_rho=%.3f | plsD1_score_rho=%.3f | plsD1_load_rho=%.3f | %.0fs",
                  density_rho,
                  run_result$pePC1_var_explained,
                  pePC1_score_rho,
                  plsD1_score_rho,
                  plsD1_loading_rho,
                  run_result$elapsed_sec))
}

# ===========================================================================
# 7. Aggregate results
# ===========================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("Aggregating results...")
message(paste(rep("=", 70), collapse = ""))

# Main summary data.frame (one row per run)
sensitivity_summary <- do.call(rbind, results_list)
rownames(sensitivity_summary) <- NULL

# Per-parameter summary (mean and range of correlations for non-default values)
parameter_summary <- sensitivity_summary |>
  dplyr::filter(status == "success") |>
  dplyr::group_by(parameter) |>
  dplyr::summarise(
    n_runs              = dplyr::n(),
    n_success           = sum(status == "success"),
    density_rho_min     = min(density_rho, na.rm = TRUE),
    density_rho_median  = median(density_rho, na.rm = TRUE),
    pePC1_var_min       = min(pePC1_var_explained, na.rm = TRUE),
    pePC1_var_max       = max(pePC1_var_explained, na.rm = TRUE),
    pePC1_rho_min       = min(pePC1_score_rho, na.rm = TRUE),
    pePC1_rho_median    = median(pePC1_score_rho, na.rm = TRUE),
    plsD1_score_rho_min = min(plsD1_score_rho, na.rm = TRUE),
    plsD1_load_rho_min  = min(plsD1_loading_rho, na.rm = TRUE),
    time_median_sec     = median(elapsed_sec, na.rm = TRUE),
    .groups = "drop"
  )

# ===========================================================================
# 8. Save outputs
# ===========================================================================

output <- list(
  sensitivity_summary = sensitivity_summary,
  parameter_summary   = parameter_summary,
  defaults            = defaults,
  param_grid          = param_grid,
  run_table           = run_table,
  session_info        = sessionInfo()
)

out_file <- file.path(out_dir, "sensitivity_results.rds")
saveRDS(object = output, file = out_file)

message(sprintf("\nResults saved to: %s", out_file))
message(sprintf("Successful runs: %d / %d",
                sum(sensitivity_summary$status == "success"),
                nrow(sensitivity_summary)))
message("\nDone.")
