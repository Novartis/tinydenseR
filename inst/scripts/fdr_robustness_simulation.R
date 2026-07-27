#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

# ──────────────────────────────────────────────────────────────────────────────
# FDR Robustness Simulation for tinydenseR
# ──────────────────────────────────────────────────────────────────────────────
#
# PURPOSE
# Evaluates empirical FDR calibration and power of the tinydenseR
# PCA-weighted q-value / density-weighted BH framework under:
#   (a) Small sample sizes  (n = 3, 4, 6, 10 per group)
#   (b) Unbalanced designs  (1:1, 2:1, 3:1 group ratios)
#   (c) Varying batch-effect magnitudes (1×, 2×, 3× default)
#
# SIMULATION DESIGN
#   Signal scenarios : 5% target cell abundance in Baseline, 2.5% in
#                      Depletion group (two-fold DA; moderate difficulty).
#   Null scenarios   : identical 5% target in both groups (no true DA).
#   Replicates       : 50 per scenario (adjustable via N_REPLICATES).
#   Total scenarios  : 4 sizes × 3 balances × 3 batch mags = 36
#                      × 2 (signal + null) = 72 scenario types.
#
# DESIGN CHOICES
#   • A local simulation function (simulate_DA_custom) mirrors the logic
#     of tinydenseR::simulate_DA_data() but adds two parameters:
#       – batch_multiplier: scales the meanlog of the log-normal batch
#         effects (default meanlog values are 0.1–0.5).
#       – null_signal: when TRUE both groups receive the same target
#         proportion, producing a true-null dataset.
#     The package source code is NOT modified.
#   • Unbalanced designs: simulate the LARGER group size for both groups,
#     then subsample the smaller group's samples (stratified by batch) when
#     loading data into the cytoset.
#   • Ground truth: each landmark's cell type is recovered by parsing
#     TDR landmark rownames (format "SampleName_event_N"), linking back
#     to cell_meta via sample name + event index.  A landmark is a
#     TRUE POSITIVE if (a) q < alpha AND (b) its originating cell belongs
#     to the "target" population.
#   • mean_cells is set to 10 000 (not the default 50 000) to keep
#     the total runtime tractable across 3 600 pipeline runs.
#     Increase MEAN_CELLS for higher-fidelity results.
#
# OUTPUT
#   A data.frame saved as an RDS file with columns:
#     scenario_id, n_per_group, balance_ratio, batch_magnitude,
#     signal_present, replicate, n_landmarks, n_significant,
#     n_TP, n_FP, empirical_FDR, power
#
# RUNTIME NOTE
#   Each RunTDR call uses .n.threads = 1 (single-threaded internally).
#   Replicates within each scenario are parallelised via mclapply (fork).
#   N_CORES workers run simultaneously; set via RhpcBLASctl::blas_get_num_procs().
#   Expect ~1-2 hours for the full grid with 8+ cores; set N_REPLICATES
#   lower (e.g. 5) for a quick check.
# ──────────────────────────────────────────────────────────────────────────────

library(tinydenseR)
library(flowCore)
library(flowWorkspace)
library(withr)
library(parallel)
library(RhpcBLASctl)

# ── Tunables ─────────────────────────────────────────────────────────────────
N_REPLICATES   <- 1L

# Detect available cores: SLURM allocation → RhpcBLASctl → detectCores
N_CORES <- max(1L, {
  # Prefer scheduler-allocated core count (most reliable on HPC)
  slurm <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)))
  lsf   <- suppressWarnings(
    as.integer(Sys.getenv("LSB_MAX_NUM_PROCESSORS", unset = NA)))
  sge   <- suppressWarnings(
    as.integer(Sys.getenv("NSLOTS", unset = NA)))

  alloc <- na.omit(c(slurm, lsf, sge))

  n <- if (length(alloc) > 0L) {
    alloc[1L]
  } else {
    # Fall back to hardware probes (max of all methods)
    max(RhpcBLASctl::blas_get_num_procs(),
        RhpcBLASctl::omp_get_num_procs(),
        RhpcBLASctl::omp_get_max_threads(),
        parallel::detectCores(logical = FALSE),
        na.rm = TRUE)
  }

  n - 1L  # reserve 1 core for the parent process
})

MEAN_CELLS     <- 10000L
SD_CELLS       <- 500L
SETTING        <- 0.05          # 5 % target-cell proportion
ALPHA          <- 0.10          # FDR threshold for calling discoveries
CL_RESOLUTION  <- 0.5
MARKERS        <- paste0("Marker", 1:5)
SEED_OFFSET    <- 1000L         # base seed; each replicate = SEED_OFFSET + rep
GROUPS         <- c("Baseline", "Depletion")
BATCHES        <- c("Batch1", "Batch2")

# ── Output directory ─────────────────────────────────────────────────────────
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

# create sub-folder for results
rd <- results_dir <-
  getwd() |>
  dirname() |>
  file.path("res",
            gsub(pattern = "\\.R$",
                 replacement = "",
                 x = script.path,
                 fixed = FALSE) |>
              basename())

if(!dir.exists(paths = rd)) dir.create(path = rd)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOCAL SIMULATION FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

#' Generate DA flow-cytometry data with configurable batch effects and null mode
#'
#' Mirrors tinydenseR::simulate_DA_data() logic with two additions:
#'   batch_multiplier – scales the meanlog of multiplicative batch shifts
#'   null_signal      – if TRUE, both groups use the same target proportion
#'
#' @return list(sample_meta, cell_meta)  identical schema to simulate_DA_data()
simulate_DA_custom <- function(groups        = GROUPS,
                               batches       = BATCHES,
                               setting       = SETTING,
                               samples_per_group = 6L,
                               mean_cells    = MEAN_CELLS,
                               sd_cells      = SD_CELLS,
                               seed          = 42L,
                               output_dir    = file.path(tempdir(), "sim_fdr"),
                               batch_multiplier = 1,
                               null_signal   = FALSE) {
  
  withr::local_preserve_seed()
  set.seed(seed)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Default batch-effect meanlog values (per marker, applied to Batch2 only)
  base_meanlog <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  base_sdlog   <- c(0.5, 0.3, 0.4, 0.3, 0.5)
  # Scale the systematic shift; keep dispersion fixed
  eff_meanlog  <- base_meanlog * batch_multiplier
  
  data_list <- list()
  
  for (group in groups) {
    for (sample_id in seq_len(samples_per_group)) {
      
      batch <- if (sample_id %% 2 == 0) "Batch2" else "Batch1"
      
      total_cells <- max(round(stats::rnorm(1, mean = mean_cells, sd = sd_cells)),
                         1000L)
      
      proportion <-
        if (null_signal) {
          setting                                       # same in both groups
        } else {
          if (group == groups[1]) setting else setting / 2
        }
      
      num_interest <- round(total_cells * proportion)
      num_other    <- total_cells - num_interest
      cell_types   <- sample(c(rep("target", num_interest),
                               rep("other",  num_other)))
      
      # Marker expression ────────────────────────────────────────────────
      marker1 <- numeric(total_cells)
      marker2 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 1.5)
      marker3 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 2.5)
      marker4 <- numeric(total_cells)
      marker5 <- numeric(total_cells)
      
      is_other  <- cell_types == "other"
      is_target <- cell_types == "target"
      
      marker1[is_other]  <- stats::rlnorm(sum(is_other),  meanlog = 0,     sdlog = 2)
      marker1[is_target] <- stats::rlnorm(sum(is_target), meanlog = 0 + 5, sdlog = 2)
      marker4[is_other]  <- stats::rlnorm(sum(is_other),  meanlog = 0,     sdlog = 1.2)
      marker4[is_target] <- stats::rlnorm(sum(is_target), meanlog = 0 + 3, sdlog = 1.2)
      marker5[is_other]  <- stats::rlnorm(sum(is_other),  meanlog = 0,     sdlog = 1.8)
      marker5[is_target] <- stats::rlnorm(sum(is_target), meanlog = 0 + 7, sdlog = 1.8)
      
      # Batch effect (Batch2 only) ──────────────────────────────────────
      if (batch == "Batch2") {
        marker1 <- marker1 * stats::rlnorm(total_cells, meanlog = eff_meanlog[1], sdlog = base_sdlog[1])
        marker2 <- marker2 * stats::rlnorm(total_cells, meanlog = eff_meanlog[2], sdlog = base_sdlog[2])
        marker3 <- marker3 * stats::rlnorm(total_cells, meanlog = eff_meanlog[3], sdlog = base_sdlog[3])
        marker4 <- marker4 * stats::rlnorm(total_cells, meanlog = eff_meanlog[4], sdlog = base_sdlog[4])
        marker5 <- marker5 * stats::rlnorm(total_cells, meanlog = eff_meanlog[5], sdlog = base_sdlog[5])
      }
      
      sample_name <- paste0(group, "_S", sample_id, "_Set", setting * 100)
      
      data_list[[length(data_list) + 1L]] <- data.frame(
        Sample    = sample_name,
        Treatment = group,
        Batch     = batch,
        Setting   = paste0(setting * 100, "%"),
        CellType  = cell_types,
        Marker1   = marker1,
        Marker2   = marker2,
        Marker3   = marker3,
        Marker4   = marker4,
        Marker5   = marker5,
        stringsAsFactors = FALSE
      )
    }
  }
  
  final_data <- do.call(rbind, data_list)
  
  # Write FCS files and build sample_meta ──────────────────────────────
  sample_names     <- unique(final_data$Sample)
  sample_meta_list <- list()
  
  for (sn in sample_names) {
    rows <- final_data[final_data$Sample == sn, ]
    mat  <- as.matrix(rows[, MARKERS])
    mat  <- log(mat)
    rownames(mat) <- paste0("event_", seq_len(nrow(mat)))
    ff   <- flowCore::flowFrame(exprs = mat)
    fcs_path <- file.path(output_dir, paste0(sn, ".fcs"))
    flowCore::write.FCS(ff, filename = fcs_path)
    
    sample_meta_list[[length(sample_meta_list) + 1L]] <- data.frame(
      Sample    = sn,
      Treatment = rows$Treatment[1],
      Batch     = rows$Batch[1],
      Setting   = rows$Setting[1],
      fcs_path  = fcs_path,
      stringsAsFactors = FALSE
    )
  }
  
  sample_meta <- do.call(rbind, sample_meta_list)
  cell_meta   <- final_data[, c("Sample", "Treatment", "Batch",
                                "Setting", "CellType")]
  
  list(sample_meta = sample_meta, cell_meta = cell_meta)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. GROUND-TRUTH MAPPING
# ══════════════════════════════════════════════════════════════════════════════

#' Map each landmark to its ground-truth cell type
#'
#' Landmark rownames have the form  "SampleName_event_N"  where N is the
#' 1-based row index within the sample's expression matrix.  The cell_meta
#' rows are in the same order as the FCS events within each sample.
#'
#' @param tdr   A TDRObj after RunTDR
#' @param cell_meta  data.frame returned by simulate_DA_custom()$cell_meta
#' @param sample_name_map  Named character vector mapping cytoset sample names
#'        (possibly with .fcs) to cell_meta$Sample names (without .fcs).
#' @return Character vector (length = n_landmarks) of "target" or "other".
map_landmark_ground_truth <- function(tdr, cell_meta, sample_name_map) {
  
  lm_names <- rownames(tdr@assay$expr)
  n_lm     <- length(lm_names)
  gt       <- character(n_lm)
  
  # Pre-split cell_meta by sample for fast lookup
  cm_by_sample <- split(cell_meta, cell_meta$Sample)
  
  # Pre-compute per-sample row indices (1-based event order)
  for (samp in names(cm_by_sample)) {
    cm_by_sample[[samp]]$.event_idx <- seq_len(nrow(cm_by_sample[[samp]]))
  }
  
  for (i in seq_len(n_lm)) {
    nm <- lm_names[i]
    
    # Extract event number: last occurrence of "_event_" followed by digits
    event_match <- regmatches(nm, regexpr("_event_\\d+$", nm))
    
    if (length(event_match) == 1L && nchar(event_match) > 0L) {
      event_idx   <- as.integer(sub("_event_", "", event_match))
      cs_sample   <- sub("_event_\\d+$", "", nm)
    } else {
      # Fallback: last underscore + digits
      event_idx <- as.integer(sub(".*_(\\d+)$", "\\1", nm))
      cs_sample <- sub("_\\d+$", "", nm)
    }
    
    # Resolve cytoset sample name → cell_meta sample name
    meta_sample <- if (cs_sample %in% names(sample_name_map)) {
      sample_name_map[[cs_sample]]
    } else {
      cs_sample
    }
    
    srows <- cm_by_sample[[meta_sample]]
    if (is.null(srows) || event_idx < 1L || event_idx > nrow(srows)) {
      gt[i] <- NA_character_
    } else {
      gt[i] <- srows$CellType[event_idx]
    }
  }
  
  gt
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. SINGLE-REPLICATE PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

#' Run one replicate: simulate → TDR → extract metrics
#'
#' @param n_group1  Integer sample count for group 1 (Baseline)
#' @param n_group2  Integer sample count for group 2 (Depletion)
#' @param batch_multiplier  Numeric batch-effect scaling factor
#' @param null_signal  Logical; TRUE for null (no DA) scenario
#' @param rep_seed  Integer seed for this replicate
#' @param alpha     Numeric FDR threshold
#'
#' @return A one-row data.frame with metric columns, or NULL on error.
run_one_replicate <- function(n_group1,
                              n_group2,
                              batch_multiplier,
                              null_signal,
                              rep_seed,
                              alpha = ALPHA) {
  
  # Temporary output directory for FCS files
  tmp_dir <- file.path(tempdir(), paste0("fdr_sim_", rep_seed))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  
  # ── Simulate data ──────────────────────────────────────────────────
  n_max <- max(n_group1, n_group2)
  
  sim <- tryCatch(
    simulate_DA_custom(
      samples_per_group = n_max,
      seed              = rep_seed,
      output_dir        = tmp_dir,
      batch_multiplier  = batch_multiplier,
      null_signal       = null_signal
    ),
    error = function(e) { message("  sim error: ", conditionMessage(e)); NULL }
  )
  if (is.null(sim)) return(NULL)
  
  # ── Subsample for imbalance ────────────────────────────────────────
  # Keep all n_group1 samples from group 1, subsample group 2 to n_group2
  meta <- sim$sample_meta
  
  keep_g1 <- meta$Sample[meta$Treatment == GROUPS[1]]
  all_g2  <- meta$Sample[meta$Treatment == GROUPS[2]]
  
  if (n_group2 < length(all_g2)) {
    # Stratified subsample by batch to minimise confounding
    g2_batches <- meta$Batch[match(all_g2, meta$Sample)]
    # Try to keep proportional batch representation
    withr::local_seed(rep_seed + 99999L)
    keep_g2 <- character(0)
    for (b in unique(g2_batches)) {
      in_batch <- all_g2[g2_batches == b]
      n_take   <- max(1L, round(n_group2 * sum(g2_batches == b) / length(all_g2)))
      n_take   <- min(n_take, length(in_batch))
      keep_g2  <- c(keep_g2, sample(in_batch, n_take))
    }
    # Trim or pad to exactly n_group2
    if (length(keep_g2) > n_group2) {
      keep_g2 <- keep_g2[seq_len(n_group2)]
    } else if (length(keep_g2) < n_group2) {
      remaining <- setdiff(all_g2, keep_g2)
      keep_g2   <- c(keep_g2, sample(remaining,
                                     min(n_group2 - length(keep_g2),
                                         length(remaining))))
    }
  } else {
    keep_g2 <- all_g2[seq_len(n_group2)]
  }
  
  # Similarly subsample group 1 if needed
  if (n_group1 < length(keep_g1)) {
    g1_batches <- meta$Batch[match(keep_g1, meta$Sample)]
    withr::local_seed(rep_seed + 88888L)
    tmp_g1 <- character(0)
    for (b in unique(g1_batches)) {
      in_batch <- keep_g1[g1_batches == b]
      n_take   <- max(1L, round(n_group1 * sum(g1_batches == b) / length(keep_g1)))
      n_take   <- min(n_take, length(in_batch))
      tmp_g1   <- c(tmp_g1, sample(in_batch, n_take))
    }
    if (length(tmp_g1) > n_group1) {
      tmp_g1 <- tmp_g1[seq_len(n_group1)]
    } else if (length(tmp_g1) < n_group1) {
      remaining <- setdiff(keep_g1, tmp_g1)
      tmp_g1 <- c(tmp_g1, sample(remaining,
                                 min(n_group1 - length(tmp_g1),
                                     length(remaining))))
    }
    keep_g1 <- tmp_g1
  }
  
  keep_samples <- c(keep_g1, keep_g2)
  meta_sub     <- meta[meta$Sample %in% keep_samples, ]
  
  # Check minimum design requirements
  n_total <- nrow(meta_sub)
  n_batches_present <- length(unique(meta_sub$Batch))
  n_groups_present  <- length(unique(meta_sub$Treatment))
  
  if (n_groups_present < 2L || n_total < 3L) {
    message("  skip: degenerate design (n=", n_total, ")")
    return(NULL)
  }
  
  # Filter cell_meta to kept samples
  cell_meta_sub <- sim$cell_meta[sim$cell_meta$Sample %in% keep_samples, ]
  
  # ── Load into cytoset ──────────────────────────────────────────────
  cs <- tryCatch({
    files <- stats::setNames(meta_sub$fcs_path, meta_sub$Sample)
    flowWorkspace::load_cytoset_from_fcs(files = files)
  }, error = function(e) { message("  cs error: ", conditionMessage(e)); NULL })
  if (is.null(cs)) return(NULL)
  
  # Resolve actual sample names in the cytoset (may include .fcs)
  cs_names <- flowWorkspace::sampleNames(cs)
  
  # Build sample-name mapping: cs_name → cell_meta$Sample
  sample_name_map <- stats::setNames(
    meta_sub$Sample,
    meta_sub$Sample  # tentative: assume names match
  )
  # If cytoset appended .fcs, update the mapping
  if (any(grepl("\\.fcs$", cs_names))) {
    sample_name_map <- stats::setNames(
      meta_sub$Sample,
      paste0(meta_sub$Sample, ".fcs")
    )
  }
  
  # Set pData
  flowWorkspace::pData(cs)$Sample    <- cs_names
  flowWorkspace::pData(cs)$Treatment <- meta_sub$Treatment[
    match(cs_names,
          if (any(grepl("\\.fcs$", cs_names)))
            paste0(meta_sub$Sample, ".fcs")
          else
            meta_sub$Sample)]
  flowWorkspace::pData(cs)$Batch <- meta_sub$Batch[
    match(cs_names,
          if (any(grepl("\\.fcs$", cs_names)))
            paste0(meta_sub$Sample, ".fcs")
          else
            meta_sub$Sample)]
  
  # ── RunTDR ─────────────────────────────────────────────────────────
  tdr <- tryCatch({
    tinydenseR::RunTDR(
      x               = cs,
      .sample.var     = "Sample",
      .harmony.var    = if (n_batches_present > 1L) "Batch" else NULL,
      .assay.type     = "cyto",
      .markers        = MARKERS,
      .seed           = rep_seed,
      .verbose        = FALSE,
      .n.threads      = 1,
      .cl.resolution.parameter = CL_RESOLUTION
    )
  }, error = function(e) { message("  TDR error: ", conditionMessage(e)); NULL })
  if (is.null(tdr)) return(NULL)
  
  # ── Linear model ───────────────────────────────────────────────────
  tdr_meta <- tdr@metadata
  
  # Build design: include Batch only if >1 batch present
  if (n_batches_present > 1L) {
    design <- model.matrix(~ Treatment + Batch, data = tdr_meta)
  } else {
    design <- model.matrix(~ Treatment, data = tdr_meta)
  }
  colnames(design) <- gsub("^Treatment|^Batch", "", colnames(design))
  
  # Identify the coefficient for the non-baseline group
  coef_name <- gsub("^Treatment", "", GROUPS[2])  # "Depletion"
  
  tdr <- tryCatch({
    tinydenseR::get.lm(x = tdr, .design = design, .verbose = FALSE)
  }, error = function(e) { message("  lm error: ", conditionMessage(e)); NULL })
  if (is.null(tdr)) return(NULL)
  
  # ── Extract logFC direction ────────────────────────────────────────
  logFC <- tryCatch(
    tdr@results$lm[["default"]]$fit$coefficients[, coef_name],
    error = function(e) NULL
  )
  if (is.null(logFC)) {
    message("  logFC extraction failed")
    return(NULL)
  }

  # ── Extract q-values ───────────────────────────────────────────────
  q_pca <- tryCatch(
    tdr@results$lm[["default"]]$fit$pca.weighted.q[, coef_name],
    error = function(e) NULL
  )
  if (is.null(q_pca)) {
    message("  q-value extraction failed")
    return(NULL)
  }
  
  # ── Ground truth ───────────────────────────────────────────────────
  gt <- map_landmark_ground_truth(tdr, cell_meta_sub, sample_name_map)
  
  is_sig    <- !is.na(q_pca) & q_pca < alpha
  is_target <- gt == "target"
  is_neg    <- !is.na(logFC) & logFC < 0
  
  n_significant <- sum(is_sig, na.rm = TRUE)
  n_TP <- sum(is_sig & is_target, na.rm = TRUE)
  n_FP <- sum(is_sig & !is_target, na.rm = TRUE)
  
  # Four-way directional breakdown of significant landmarks
  n_sig_target_neg <- sum(is_sig & is_target  & is_neg,  na.rm = TRUE)
  n_sig_target_pos <- sum(is_sig & is_target  & !is_neg, na.rm = TRUE)
  n_sig_other_neg  <- sum(is_sig & !is_target & is_neg,  na.rm = TRUE)
  n_sig_other_pos  <- sum(is_sig & !is_target & !is_neg, na.rm = TRUE)
  
  # Total target landmarks (for power denominator)
  n_target_landmarks <- sum(is_target, na.rm = TRUE)
  
  empirical_FDR <- if (n_significant > 0L) n_FP / n_significant else NA_real_
  power         <- if (n_target_landmarks > 0L) n_TP / n_target_landmarks else NA_real_
  
  data.frame(
    n_landmarks      = length(q_pca),
    n_target_lm      = n_target_landmarks,
    n_significant    = n_significant,
    n_TP             = n_TP,
    n_FP             = n_FP,
    n_sig_target_neg = n_sig_target_neg,
    n_sig_target_pos = n_sig_target_pos,
    n_sig_other_neg  = n_sig_other_neg,
    n_sig_other_pos  = n_sig_other_pos,
    empirical_FDR    = empirical_FDR,
    power            = power,
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. SCENARIO GRID
# ══════════════════════════════════════════════════════════════════════════════

build_scenario_grid <- function() {
  n_per_group     <- c(3L, 4L, 6L, 10L)
  balance_ratio   <- c(1L, 2L, 3L)            # group1:group2 ratio
  batch_magnitude <- c(1, 2, 3)               # multiplier on meanlog
  signal_present  <- c(TRUE, FALSE)
  
  grid <- expand.grid(
    n_per_group     = n_per_group,
    balance_ratio   = balance_ratio,
    batch_magnitude = batch_magnitude,
    signal_present  = signal_present,
    stringsAsFactors = FALSE
  )
  
  # Compute per-group sample sizes
  # n_per_group = size of the larger group (group 1 = Baseline)
  # group 2 = max(2, round(n_per_group / balance_ratio))
  grid$n_group1 <- grid$n_per_group
  grid$n_group2 <- pmax(2L, as.integer(round(grid$n_per_group / grid$balance_ratio)))
  
  # Special case: ratio 1 means equal groups
  
  grid$n_group2[grid$balance_ratio == 1L] <-
    grid$n_group1[grid$balance_ratio == 1L]
  
  grid$scenario_id <- seq_len(nrow(grid))
  
  grid
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. MAIN SIMULATION LOOP
# ══════════════════════════════════════════════════════════════════════════════

run_all_simulations <- function(n_replicates = N_REPLICATES,
                                n_cores      = N_CORES,
                                output_dir   = results_dir) {
  
  scenarios <- build_scenario_grid()
  n_scen    <- nrow(scenarios)
  
  message("Running ", n_scen, " scenarios × ", n_replicates,
          " replicates = ", n_scen * n_replicates, " total runs")

  message("Parallel workers: ", n_cores,
          " (mclapply fork, detected via RhpcBLASctl)")
  
  all_results <- vector("list", n_scen)
  
  for (s in seq_len(n_scen)) {
    sc <- scenarios[s, ]
    
    message(sprintf(
      "\n── Scenario %d/%d: n=%d, ratio=%d:1, batch=%.0fx, signal=%s ──",
      s, n_scen, sc$n_per_group, sc$balance_ratio,
      sc$batch_magnitude, if (sc$signal_present) "YES" else "NULL"
    ))
    
    # Build replicate seeds
    rep_seeds <- SEED_OFFSET + (s - 1L) * n_replicates + seq_len(n_replicates)
    
    # Run replicates in parallel via mclapply
    rep_results <- parallel::mclapply(
      X = seq_len(n_replicates),
      FUN = function(r) {
        
        # Suppress BLAS multi-threading inside each forked worker
        # to avoid thread oversubscription (N_CORES workers already saturate CPUs)
        invisible(utils::capture.output({
          RhpcBLASctl::blas_set_num_threads(1L)
          RhpcBLASctl::omp_set_num_threads(1L)
        }, type = "output"))
        
        rep_seed <- rep_seeds[r]
        
        metrics <- tryCatch(
          run_one_replicate(
            n_group1         = sc$n_group1,
            n_group2         = sc$n_group2,
            batch_multiplier = sc$batch_magnitude,
            null_signal      = !sc$signal_present,
            rep_seed         = rep_seed,
            alpha            = ALPHA
          ),
          error = function(e) {
            message("  UNEXPECTED error (rep ", r, "): ", conditionMessage(e))
            NULL
          }
        )
        
        if (is.null(metrics)) {
          metrics <- data.frame(
            n_landmarks      = NA_integer_,
            n_target_lm      = NA_integer_,
            n_significant    = NA_integer_,
            n_TP             = NA_integer_,
            n_FP             = NA_integer_,
            n_sig_target_neg = NA_integer_,
            n_sig_target_pos = NA_integer_,
            n_sig_other_neg  = NA_integer_,
            n_sig_other_pos  = NA_integer_,
            empirical_FDR    = NA_real_,
            power            = NA_real_,
            stringsAsFactors = FALSE
          )
        }
        
        # Annotate with scenario info
        metrics$scenario_id      <- sc$scenario_id
        metrics$n_per_group      <- sc$n_per_group
        metrics$n_group1         <- sc$n_group1
        metrics$n_group2         <- sc$n_group2
        metrics$balance_ratio    <- paste0(sc$balance_ratio, ":1")
        metrics$batch_magnitude  <- sc$batch_magnitude
        metrics$signal_present   <- sc$signal_present
        metrics$replicate        <- r
        metrics$seed             <- rep_seed
        
        metrics
      },
      mc.cores = n_cores,
      mc.set.seed = FALSE  # we manage seeds explicitly per replicate
    )
    
    # Collect results for this scenario
    # mclapply returns a list; check for try-error objects from crashed workers
    rep_results <- lapply(rep_results, function(x) {
      if (inherits(x, "try-error")) {
        message("  Worker crash: ", attr(x, "condition")$message)
        data.frame(
          n_landmarks = NA_integer_, n_target_lm = NA_integer_,
          n_significant = NA_integer_, n_TP = NA_integer_,
          n_FP = NA_integer_,
          n_sig_target_neg = NA_integer_, n_sig_target_pos = NA_integer_,
          n_sig_other_neg = NA_integer_, n_sig_other_pos = NA_integer_,
          empirical_FDR = NA_real_, power = NA_real_,
          scenario_id = sc$scenario_id, n_per_group = sc$n_per_group,
          n_group1 = sc$n_group1, n_group2 = sc$n_group2,
          balance_ratio = paste0(sc$balance_ratio, ":1"),
          batch_magnitude = sc$batch_magnitude, signal_present = sc$signal_present,
          replicate = NA_integer_, seed = NA_integer_,
          stringsAsFactors = FALSE
        )
      } else {
        x
      }
    })
    
    all_results[[s]] <- do.call(rbind, rep_results)
    
    message(sprintf("  completed %d replicates (%d OK, %d failed)",
                    n_replicates,
                    sum(!is.na(all_results[[s]]$n_significant)),
                    sum(is.na(all_results[[s]]$n_significant))))
    
    # Checkpoint: save intermediate results after each scenario
    intermediate <- do.call(rbind, Filter(Negate(is.null), all_results))
    saveRDS(intermediate,
            file = file.path(output_dir, "fdr_robustness_results_checkpoint.rds"))
  }
  
  results <- do.call(rbind, all_results)
  
  # Reorder columns
  col_order <- c("scenario_id", "n_per_group", "n_group1", "n_group2",
                 "balance_ratio", "batch_magnitude", "signal_present",
                 "replicate", "seed", "n_landmarks", "n_target_lm",
                 "n_significant", "n_TP", "n_FP",
                 "n_sig_target_neg", "n_sig_target_pos",
                 "n_sig_other_neg", "n_sig_other_pos",
                 "empirical_FDR", "power")
  results <- results[, col_order]
  rownames(results) <- NULL
  
  results
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. RESULT AGGREGATION
# ══════════════════════════════════════════════════════════════════════════════

#' Summarise simulation results per scenario
#'
#' @param results  data.frame from run_all_simulations()
#' @return data.frame with one row per scenario and summary statistics
summarise_results <- function(results) {
  
  scenarios <- unique(results[, c("scenario_id", "n_per_group", "n_group1",
                                  "n_group2", "balance_ratio",
                                  "batch_magnitude", "signal_present")])
  
  do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    sc  <- scenarios[i, ]
    sub <- results[results$scenario_id == sc$scenario_id, ]
    
    # Remove NA replicates
    sub_valid <- sub[!is.na(sub$n_significant), ]
    n_valid   <- nrow(sub_valid)
    
    data.frame(
      scenario_id      = sc$scenario_id,
      n_per_group      = sc$n_per_group,
      n_group1         = sc$n_group1,
      n_group2         = sc$n_group2,
      balance_ratio    = sc$balance_ratio,
      batch_magnitude  = sc$batch_magnitude,
      signal_present   = sc$signal_present,
      n_replicates     = n_valid,
      n_failed         = nrow(sub) - n_valid,
      
      # Mean empirical FDR across replicates (excluding NAs from 0-discovery reps)
      mean_empirical_FDR = mean(sub_valid$empirical_FDR, na.rm = TRUE),
      sd_empirical_FDR   = sd(sub_valid$empirical_FDR, na.rm = TRUE),
      median_empirical_FDR = median(sub_valid$empirical_FDR, na.rm = TRUE),
      
      # Mean power (signal scenarios only)
      mean_power       = mean(sub_valid$power, na.rm = TRUE),
      sd_power         = sd(sub_valid$power, na.rm = TRUE),
      median_power     = median(sub_valid$power, na.rm = TRUE),
      
      # Discovery counts
      mean_n_sig       = mean(sub_valid$n_significant, na.rm = TRUE),
      mean_n_TP        = mean(sub_valid$n_TP, na.rm = TRUE),
      mean_n_FP        = mean(sub_valid$n_FP, na.rm = TRUE),

      # Directional breakdown (diagnostic: inspect before reclassifying)
      mean_sig_target_neg = mean(sub_valid$n_sig_target_neg, na.rm = TRUE),
      mean_sig_target_pos = mean(sub_valid$n_sig_target_pos, na.rm = TRUE),
      mean_sig_other_neg  = mean(sub_valid$n_sig_other_neg, na.rm = TRUE),
      mean_sig_other_pos  = mean(sub_valid$n_sig_other_pos, na.rm = TRUE),
      
      # Type I error rate (null scenarios): fraction of reps with ≥1 discovery
      type1_any_disc   = mean(sub_valid$n_significant > 0, na.rm = TRUE),
      
      # FDR calibration: is mean empirical FDR ≤ nominal alpha?
      fdr_calibrated   = ifelse(
        sc$signal_present,
        mean(sub_valid$empirical_FDR, na.rm = TRUE) <= ALPHA,
        mean(sub_valid$n_significant > 0, na.rm = TRUE) <= ALPHA
      ),
      
      stringsAsFactors = FALSE
    )
  }))
}

# ══════════════════════════════════════════════════════════════════════════════
# 7. MAIN ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

main <- function() {
  message("=== FDR Robustness Simulation ===")
  message("Replicates : ", N_REPLICATES)
  message("Cores      : ", N_CORES, " (mclapply fork workers)")
  message("Mean cells : ", MEAN_CELLS)
  message("Setting    : ", SETTING * 100, "% target proportion")
  message("Alpha      : ", ALPHA)
  message("Output     : ", results_dir)
  
  results <- run_all_simulations()
  
  # Save full results
  out_path <- file.path(results_dir, "fdr_robustness_results.rds")
  saveRDS(results, file = out_path)
  message("\nFull results saved to: ", out_path)
  
  # Save summary
  summary_df  <- summarise_results(results)
  summary_path <- file.path(results_dir, "fdr_robustness_summary.rds")
  saveRDS(summary_df, file = summary_path)
  message("Summary saved to: ", summary_path)
  
  # Also write summary as CSV for quick inspection
  csv_path <- file.path(results_dir, "fdr_robustness_summary.csv")
  utils::write.csv(summary_df, file = csv_path, row.names = FALSE)
  message("Summary CSV: ", csv_path)
  
  # Print summary
  message("\n── Summary ──")
  print(summary_df[, c("scenario_id", "n_per_group", "balance_ratio",
                       "batch_magnitude", "signal_present", "n_replicates",
                       "mean_empirical_FDR", "mean_power",
                       "type1_any_disc", "fdr_calibrated")])
  
  invisible(results)
}

Sys.setenv(FDR_SIM_RUN = "1")
# Run if sourced or executed directly
if (!interactive() || identical(Sys.getenv("FDR_SIM_RUN"), "1")) {
  main()
}
Sys.unsetenv("FDR_SIM_RUN")