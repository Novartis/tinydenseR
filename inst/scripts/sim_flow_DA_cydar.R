#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# Licensed under the MIT License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https://www.mit.edu/~amini/LICENSE.md
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#####

library(cydar)
library(flowCore)
library(edgeR)
library(tinydenseR)
library(tidyverse)

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
rd <-
  script.path |>
  dirname() |>
  dirname() |>
  file.path("res",
            "sim_flow_DA_cydar")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

set.seed(42)

# Simulate DA data using tinydenseR API (same seed / RNG sequence as original)
sim_data <- tinydenseR::simulate_DA_data()
final_data_DA <- sim_data$cell_meta |>
  dplyr::mutate(Treatment = factor(x = Treatment,
                                   levels = c("Baseline", "Depletion")))

settings <- c(0.005, 0.05, 0.5)

# Function to convert simulated data to cydar format
# Reads pre-written FCS files from simulate_DA_data() output.
prepare_cydar_data <- function(sample_meta, setting_filter) {

  # Filter sample metadata for specific setting
  meta_filtered <- sample_meta |>
    dplyr::filter(Setting == setting_filter)

  sample_ids <- meta_filtered$Sample

  # Read FCS files as list of matrices (cydar expects list of matrices)
  mat_list <- list()

  for (i in seq_along(sample_ids)) {
    ff <- flowCore::read.FCS(meta_filtered$fcs_path[i],
                             transformation = FALSE,
                             truncate_max_range = FALSE)
    mat_list[[i]] <- flowCore::exprs(ff)
  }

  names(mat_list) <- sample_ids

  # Create experiment_info
  experiment_info <- meta_filtered |>
    dplyr::select(Sample, Treatment, Batch) |>
    as.data.frame()

  rownames(experiment_info) <- experiment_info$Sample

  return(list(
    mat_list = mat_list,
    experiment_info = experiment_info
  ))
}

# Note: simulated FCS data is already log-transformed by simulate_DA_data().
# No additional transformation (e.g., logicle) is applied, consistent with
# the pre-normalized data format. No gating is needed as the simulation
# does not include debris, dead cells, or QC markers.

# Function to run cydar analysis
run_cydar_analysis <- function(cydar_data, setting_name) {

  cat("Running cydar analysis for", setting_name, "\n")

  # prepareCellData expects a list of numeric matrices (one per sample)
  cd <- cydar::prepareCellData(cydar_data$mat_list)

  # Count cells in hyperspheres
  # Use a tolerance appropriate for log-transformed data
  cd <- cydar::countCells(cd,
                          downsample = 10,
                          BPPARAM = BiocParallel::SerialParam())

  # Build design matrix for edgeR-based DA testing
  experiment_info <- cydar_data$experiment_info
  design <- model.matrix(~ Treatment + Batch, data = experiment_info)

  # Use edgeR for differential abundance testing on hypersphere counts
  y <- DGEList(counts = assay(cd), lib.size = cd$totals)

  # Filter low-abundance hyperspheres (cydar vignette recommendation).
  # Removes near-empty hyperspheres that inflate multiple testing burden
  # and can produce unstable dispersion estimates.
  keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
  cd <- cd[keep,]
  y <- y[keep,]

  y <- estimateDisp(y, design)

  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, coef = 2)  # Treatment effect

  # Spatial FDR correction
  qvals <- cydar::spatialFDR(intensities(cd), res$table$PValue)

  results_df <- data.frame(
    hypersphere = seq_len(nrow(res$table)),
    logFC = res$table$logFC,
    p_val = res$table$PValue,
    p_adj = qvals,
    significant = qvals < 0.1
  )

  return(list(
    cd = cd,
    results = results_df,
    fit = fit
  ))
}

# Function to assess overlap with ground truth
assess_ground_truth_overlap <- function(cydar_results, cydar_data, cell_meta,
                                        setting_filter) {

  cd <- cydar_results$cd
  results_df <- cydar_results$results

  # Get significant hyperspheres
  sig_idx <- which(results_df$significant)

  if (length(sig_idx) == 0) {
    return(list(
      n_significant = 0,
      precision = NA_real_,
      recall = NA_real_
    ))
  }

  # Extract cell assignments for significant hyperspheres
  # Build explicit cell-to-metadata mapping.
  # cydar pools cells in input order (sample 1 rows, then sample 2, etc.).
  # We reconstruct this mapping explicitly rather than relying on
  # undocumented ordering internals.
  cell_counts <- vapply(cydar_data$mat_list, nrow, integer(1))
  cell_sample_map <- rep(names(cydar_data$mat_list), cell_counts)

  setting_cells <- cell_meta |>
    dplyr::filter(Setting == setting_filter)

  # Order cells by sample in the same order as mat_list
  ordered_cells <- setting_cells |>
    dplyr::arrange(match(Sample, names(cydar_data$mat_list)))

  # Sanity check: pooled cell count must match metadata rows
  stopifnot(sum(cell_counts) == nrow(ordered_cells))

  # cellAssignments maps hyperspheres to the cells they contain
  cell_assign <- cellAssignments(cd)
  sig_cells <- unique(unlist(cell_assign[sig_idx]))

  # Check which significant cells are target
  n_total_target <- sum(ordered_cells$CellType == "target")
  n_sig_cells <- length(sig_cells)

  # Cells in significant hyperspheres that are target
  sig_cell_types <- ordered_cells$CellType[sig_cells]
  n_sig_target <- sum(sig_cell_types == "target", na.rm = TRUE)

  precision <- n_sig_target / n_sig_cells

  recall <- n_sig_target / n_total_target

  return(list(
    n_significant = sum(results_df$significant),
    n_sig_cells = n_sig_cells,
    n_sig_target = n_sig_target,
    precision = precision,
    recall = recall
  ))
}

# Function to calculate ground truth proportions for comparison
calculate_ground_truth <- function(data, setting_filter) {

  ground_truth <- data |>
    dplyr::filter(Setting == setting_filter) |>
    dplyr::group_by(Sample, Treatment, Batch) |>
    dplyr::summarise(
      total_cells = n(),
      target_cells = sum(CellType == "target"),
      target_proportion = target_cells / total_cells,
      .groups = "drop"
    )

  baseline_prop <- mean(ground_truth$target_proportion[ground_truth$Treatment == "Baseline"])
  depletion_prop <- mean(ground_truth$target_proportion[ground_truth$Treatment == "Depletion"])
  expected_logFC <- log2(depletion_prop / baseline_prop)

  cat("\nGround Truth for", setting_filter, ":\n")
  cat("Baseline mean proportion:", round(baseline_prop, 4), "\n")
  cat("Depletion mean proportion:", round(depletion_prop, 4), "\n")
  cat("Expected log2FC:", round(expected_logFC, 3), "\n\n")

  return(list(
    ground_truth = ground_truth,
    expected_logFC = expected_logFC,
    baseline_prop = baseline_prop,
    depletion_prop = depletion_prop
  ))
}

# Run analysis for each setting
results_all <- list()
ground_truth_all <- list()
timing_all <- list()

for (setting in settings) {
  setting_name <- paste0(setting * 100, "%")

  cat("\n", rep("=", 50), "\n")

  cat("Processing setting:", setting_name, "\n")
  cat(rep("=", 50), "\n")

  # Calculate ground truth
  ground_truth <- calculate_ground_truth(final_data_DA, setting_name)
  ground_truth_all[[setting_name]] <- ground_truth

  # Prepare data
  cydar_data <- prepare_cydar_data(sim_data$sample_meta, setting_name)

  # Run cydar analysis with timing
  elapsed <- system.time({
    cydar_results <- run_cydar_analysis(cydar_data, setting_name)
  })

  timing_all[[setting_name]] <- elapsed["elapsed"]

  # Assess overlap with ground truth
  overlap <- assess_ground_truth_overlap(
    cydar_results, cydar_data, final_data_DA, setting_name
  )

  # Store results
  results_all[[setting_name]] <- list(
    cydar_results = cydar_results,
    overlap = overlap,
    elapsed = elapsed["elapsed"]
  )

  cat("\ncydar results for", setting_name, ":\n")
  cat("  Significant hyperspheres (FDR < 0.1):", overlap$n_significant, "\n")
  cat("  Precision (target cells among sig cells):", round(overlap$precision, 3), "\n")
  cat("  Recall (sig target cells / total target):", round(overlap$recall, 3), "\n")
  cat("  Elapsed time (s):", round(elapsed["elapsed"], 2), "\n")
}

# Summary comparison across settings
cat("\n", rep("=", 60), "\n")
cat("SUMMARY: cydar DA detection across settings\n")
cat(rep("=", 60), "\n")

summary_comparison <- data.frame()

for (setting in settings) {
  setting_name <- paste0(setting * 100, "%")

  overlap <- results_all[[setting_name]]$overlap
  gt <- ground_truth_all[[setting_name]]

  summary_row <- data.frame(
    Setting = setting_name,
    Expected_logFC = round(gt$expected_logFC, 3),
    N_Significant_Hyperspheres = overlap$n_significant,
    Precision = round(overlap$precision, 3),
    Recall = round(overlap$recall, 3),
    Elapsed_Seconds = round(results_all[[setting_name]]$elapsed, 2),
    Baseline_Prop = round(gt$baseline_prop, 4),
    Depletion_Prop = round(gt$depletion_prop, 4)
  )

  summary_comparison <- rbind(summary_comparison, summary_row)
}

print(summary_comparison)

# Save results
saveRDS(summary_comparison, file = file.path(rd, "cydar_summary.rds"))

# =========================================================================
# Supplementary Table: cydar DA detection summary
# =========================================================================

# Compute logFC direction breakdown for significant hyperspheres.
# Due to compositionality, cydar detects both the depleted population
# (negative logFC) and the compensatorily enriched non-target population
# (positive logFC).
logfc_direction <- lapply(settings, function(setting) {
  setting_name <- paste0(setting * 100, "%")
  res <- results_all[[setting_name]]$cydar_results$results
  sig <- res[res$significant, ]
  data.frame(
    Setting = setting_name,
    n_negative = sum(sig$logFC < 0),
    n_positive = sum(sig$logFC > 0),
    median_neg_logFC = if (any(sig$logFC < 0)) median(sig$logFC[sig$logFC < 0]) else NA_real_,
    median_pos_logFC = if (any(sig$logFC > 0)) median(sig$logFC[sig$logFC > 0]) else NA_real_,
    stringsAsFactors = FALSE
  )
}) |>
  do.call(what = rbind)

cat("\nLogFC direction breakdown (compositional effect):\n")
print(logfc_direction, row.names = FALSE)

# Format for manuscript supplementary table
supp_table <- summary_comparison |>
  dplyr::left_join(logfc_direction, by = "Setting") |>
  dplyr::transmute(
    `Target frequency` = Setting,
    `Baseline (%)` = Baseline_Prop * 100,
    `Depletion (%)` = Depletion_Prop * 100,
    `Expected log2FC` = Expected_logFC,
    `Significant hyperspheres (q < 0.1)` = N_Significant_Hyperspheres,
    `Negative logFC` = n_negative,
    `Positive logFC` = n_positive,
    `Cell-level precision` = ifelse(is.na(Precision), "--", sprintf("%.1f%%", Precision * 100)),
    `Cell-level recall` = ifelse(is.na(Recall), "--", sprintf("%.1f%%", Recall * 100)),
    `Time (s)` = Elapsed_Seconds
  )

cat("\n", rep("=", 60), "\n")
cat("TABLE: cydar summary for manuscript\n")
cat(rep("=", 60), "\n")
print(supp_table, row.names = FALSE)

write.csv(supp_table, file = file.path(rd, "cydar_supp_table.csv"),
          row.names = FALSE)

# =========================================================================
# Diagnostic bar plot: significant hypersphere logFC distributions
# =========================================================================

logfc_data <- lapply(settings, function(setting) {
  setting_name <- paste0(setting * 100, "%")
  res <- results_all[[setting_name]]$cydar_results$results
  if (sum(res$significant) == 0) return(NULL)
  data.frame(
    Setting = setting_name,
    logFC = res$logFC[res$significant],
    stringsAsFactors = FALSE
  )
}) |>
  (\(x) do.call(rbind, x[!vapply(x, is.null, logical(1))]))()

if (!is.null(logfc_data) && nrow(logfc_data) > 0) {
  p <- ggplot2::ggplot(logfc_data,
                       ggplot2::aes(x = logFC)) +
    ggplot2::geom_histogram(bins = 50, fill = "grey30", color = "white") +
    ggplot2::facet_wrap(~ Setting, scales = "free_y") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(x = "log2 fold change (Depletion vs Baseline)",
                  y = "Significant hyperspheres",
                  title = "cydar: logFC distribution of significant hyperspheres") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(plot = p,
                  filename = file.path(rd, "cydar_sig_logFC_histogram.png"),
                  width = 6, height = 3, dpi = 300, bg = "white")
  rm(p)
}

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat("cydar hypersphere-based DA testing on simulated data\n")
cat("Results saved to:", rd, "\n")
cat(rep("=", 60), "\n")

sessionInfo()
