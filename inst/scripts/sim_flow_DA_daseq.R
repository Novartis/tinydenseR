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

library(DAseq)
library(flowCore)
library(tinydenseR)
library(tidyverse)
library(ggpubr)

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
rd <- file.path(dirname(script.path), "results", "sim_flow_DA_daseq")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

set.seed(42)

# Simulate DA data using tinydenseR API (same seed / RNG sequence as original)
sim_data <- tinydenseR::simulate_DA_data(mean_cells = 5000)
final_data_DA <- sim_data$cell_meta |>
  dplyr::mutate(Treatment = factor(x = Treatment,
                                   levels = c("Baseline", "Depletion")))

settings <- c(0.005, 0.05, 0.5)

# Function to prepare data for DA-seq
# DA-seq operates on a shared embedding with per-cell sample labels.
# cell.labels must be sample IDs; labels.1 / labels.2 are the sample names
# belonging to each condition (see CompCy-lab/benchmarkDA).
prepare_daseq_data <- function(sample_meta, cell_meta, setting_filter) {
  
  # Filter to setting
  meta_filtered <- sample_meta |>
    dplyr::filter(Setting == setting_filter)
  
  cells_filtered <- cell_meta |>
    dplyr::filter(Setting == setting_filter)
  
  # Read FCS files and concatenate expression matrices
  expr_list <- list()
  
  for (i in seq_len(nrow(meta_filtered))) {
    ff <- flowCore::read.FCS(meta_filtered$fcs_path[i],
                             transformation = FALSE,
                             truncate_max_range = FALSE)
    expr_list[[i]] <- flowCore::exprs(ff)
  }
  
  expr_mat <- do.call(rbind, expr_list)
  
  # DA-seq needs per-cell sample IDs as cell.labels
  sample_labels <- as.character(cells_filtered$Sample)
  
  # Sample names belonging to each condition
  samples_baseline <- meta_filtered$Sample[meta_filtered$Treatment == "Baseline"]
  samples_depletion <- meta_filtered$Sample[meta_filtered$Treatment == "Depletion"]
  
  return(list(
    embedding = expr_mat,
    sample_labels = sample_labels,
    samples_1 = samples_baseline,
    samples_2 = samples_depletion,
    cell_meta = cells_filtered
  ))
}

# Function to run DA-seq analysis
# Following CompCy-lab/benchmarkDA: use getDAcells with sample-level labels.
# getDAcells internally computes permutation-based thresholds and returns
# da.up / da.down with those thresholds already applied.
# The threshold values themselves are NOT stored as $pred.thres — they must
# be recovered from $rand.pred (the null predictions from permuted labels).
run_daseq_analysis <- function(daseq_data, setting_name) {
  
  cat("Running DA-seq analysis for", setting_name, "\n")
  
  embedding <- daseq_data$embedding
  
  # DA-seq: score cells by local enrichment across k values.
  # cell.labels = per-cell sample IDs
  # labels.1 / labels.2 = sample names per condition
  # Internally, getDAcells:
  #   1. Computes multi-scale k-NN ratio scores
  #   2. Trains logistic regression on scores
  #   3. Runs permutations to get null score distribution
  #   4. Sets threshold = [min(null), max(null)]
  #   5. Returns da.up / da.down with correct thresholds applied
  da_cells <- getDAcells(
    X = embedding,
    cell.labels = daseq_data$sample_labels,
    labels.1 = daseq_data$samples_1,
    labels.2 = daseq_data$samples_2,
    plot.embedding = embedding[, 1:2],
    size = 1,
    do.plot = FALSE
  )
  
  # Extract the permutation-derived thresholds from rand.pred
  # (getDAcells does NOT return $pred.thres — it stores null predictions
  # in $rand.pred and applies [min, max] as thresholds internally)
  rand_preds <- unlist(da_cells$rand.pred)
  pred_thres <- c(min(rand_preds), max(rand_preds))
  cat("  Permutation thresholds: [", round(pred_thres[1], 4), ",",
      round(pred_thres[2], 4), "]\n")
  cat("  DA cells from getDAcells — up:", length(da_cells$da.up),
      " down:", length(da_cells$da.down), "\n")
  
  return(list(
    da_cells = da_cells,
    pred_thres = pred_thres
  ))
}

# Function to assess overlap with ground truth
assess_ground_truth_overlap <- function(daseq_results, daseq_data,
                                        setting_name) {
  
  da_cells <- daseq_results$da_cells
  cell_meta <- daseq_data$cell_meta
  
  # da.up = cells enriched in condition 2 (Depletion)
  # da.down = cells enriched in condition 1 (Baseline), i.e. depleted in cond 2
  # These are set by getDAcells using permutation-derived thresholds
  # Since target cells are depleted in Depletion, they should appear in da.down
  da_cell_idx_up <- da_cells$da.up
  da_cell_idx_down <- da_cells$da.down
  sig_cells <- da_cell_idx_down
  
  n_sig <- length(sig_cells)
  
  if (n_sig == 0) {
    return(list(
      n_da_cells = 0,
      n_da_cells_up = length(da_cell_idx_up),
      n_da_cells_down = 0,
      precision = NA_real_,
      recall = NA_real_
    ))
  }
  
  # Check overlap with ground truth target cells
  cell_types <- cell_meta$CellType
  n_total_target <- sum(cell_types == "target")
  
  sig_cell_types <- cell_types[sig_cells]
  n_sig_target <- sum(sig_cell_types == "target")
  
  precision <- n_sig_target / n_sig
  recall <- n_sig_target / n_total_target
  
  return(list(
    n_da_cells = n_sig,
    n_da_cells_up = length(da_cell_idx_up),
    n_da_cells_down = length(da_cell_idx_down),
    n_sig_target = n_sig_target,
    precision = precision,
    recall = recall
  ))
}

# Function to calculate ground truth proportions
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

for (setting in settings) {
  setting_name <- paste0(setting * 100, "%")
  
  cat("\n", rep("=", 50), "\n")
  cat("Processing setting:", setting_name, "\n")
  cat(rep("=", 50), "\n")
  
  # Calculate ground truth
  ground_truth <- calculate_ground_truth(final_data_DA, setting_name)
  ground_truth_all[[setting_name]] <- ground_truth
  
  # Prepare data
  daseq_data <- prepare_daseq_data(sim_data$sample_meta,
                                   final_data_DA,
                                   setting_name)
  
  # Run DA-seq analysis with timing
  elapsed <- system.time({
    daseq_results <- run_daseq_analysis(daseq_data, setting_name)
  })
  
  # Assess overlap with ground truth
  overlap <- assess_ground_truth_overlap(daseq_results, daseq_data,
                                         setting_name)
  
  # Store results
  results_all[[setting_name]] <- list(
    daseq_results = daseq_results,
    overlap = overlap,
    elapsed = elapsed["elapsed"]
  )
  
  cat("\nDA-seq results for", setting_name, ":\n")
  cat("  DA cells (depleted in Depletion):", overlap$n_da_cells, "\n")
  cat("  DA cells (enriched in Depletion):", overlap$n_da_cells_up, "\n")
  cat("  Precision (target cells among DA cells):", round(overlap$precision, 3), "\n")
  cat("  Recall (DA target cells / total target):", round(overlap$recall, 3), "\n")
  cat("  Elapsed time (s):", round(elapsed["elapsed"], 2), "\n")
}

# Summary comparison across settings
cat("\n", rep("=", 60), "\n")
cat("SUMMARY: DA-seq DA detection across settings\n")
cat(rep("=", 60), "\n")

summary_comparison <- data.frame()

for (setting in settings) {
  setting_name <- paste0(setting * 100, "%")
  
  overlap <- results_all[[setting_name]]$overlap
  gt <- ground_truth_all[[setting_name]]
  
  summary_row <- data.frame(
    Setting = setting_name,
    Expected_logFC = round(gt$expected_logFC, 3),
    N_DA_Cells_Down = overlap$n_da_cells,
    N_DA_Cells_Up = overlap$n_da_cells_up,
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
saveRDS(summary_comparison, file = file.path(rd, "daseq_summary.rds"))

# =========================================================================
# Supplementary Table: DA-seq DA detection summary
# =========================================================================

# Compute DA score breakdown for each setting.
# DA-seq assigns per-cell scores; cells exceeding the permutation-derived
# threshold are classified as DA (up = enriched in Depletion, down = depleted).
score_summary <- lapply(settings, function(setting) {
  setting_name <- paste0(setting * 100, "%")
  res <- results_all[[setting_name]]$daseq_results
  da_pred <- res$da_cells$da.pred
  thres <- res$pred_thres
  data.frame(
    Setting = setting_name,
    n_cells = length(da_pred),
    median_score = round(median(da_pred), 4),
    mean_score = round(mean(da_pred), 4),
    threshold_low = round(thres[1], 4),
    threshold_high = round(thres[2], 4),
    stringsAsFactors = FALSE
  )
}) |>
  do.call(what = rbind)

cat("\nDA score summary per setting:\n")
print(score_summary, row.names = FALSE)

# Format for manuscript supplementary table
supp_table <- summary_comparison |>
  dplyr::left_join(score_summary, by = "Setting") |>
  dplyr::transmute(
    `Target frequency` = Setting,
    `Baseline (%)` = Baseline_Prop * 100,
    `Depletion (%)` = Depletion_Prop * 100,
    `Expected log2FC` = Expected_logFC,
    `DA cells (depleted)` = N_DA_Cells_Down,
    `DA cells (enriched)` = N_DA_Cells_Up,
    `Score threshold` = sprintf("[%.3f, %.3f]", threshold_low, threshold_high),
    `Cell-level precision` = ifelse(is.na(Precision), "--", sprintf("%.1f%%", Precision * 100)),
    `Cell-level recall` = ifelse(is.na(Recall), "--", sprintf("%.1f%%", Recall * 100)),
    `Time (s)` = Elapsed_Seconds
  )

cat("\n", rep("=", 60), "\n")
cat("TABLE: DA-seq summary for manuscript\n")
cat(rep("=", 60), "\n")
print(supp_table, row.names = FALSE)

write.csv(supp_table, file = file.path(rd, "daseq_supp_table.csv"),
          row.names = FALSE)

# =========================================================================
# Diagnostic plot: DA score distributions per setting
# =========================================================================

score_data <- lapply(settings, function(setting) {
  setting_name <- paste0(setting * 100, "%")
  da_pred <- results_all[[setting_name]]$daseq_results$da_cells$da.pred
  data.frame(
    Setting = setting_name,
    da_score = da_pred,
    stringsAsFactors = FALSE
  )
}) |>
  do.call(what = rbind)

if (nrow(score_data) > 0) {
  # Add threshold lines per setting
  thres_df <- do.call(rbind, lapply(settings, function(setting) {
    setting_name <- paste0(setting * 100, "%")
    thres <- results_all[[setting_name]]$daseq_results$pred_thres
    data.frame(
      Setting = setting_name,
      xintercept = thres,
      stringsAsFactors = FALSE
    )
  }))
  
  p <- ggplot2::ggplot(score_data,
                       ggplot2::aes(x = da_score)) +
    ggplot2::geom_histogram(bins = 80, fill = "grey30", color = "white") +
    ggplot2::facet_wrap(~ Setting, scales = "free_y") +
    ggplot2::geom_vline(data = thres_df,
                        ggplot2::aes(xintercept = xintercept),
                        linetype = "dashed", color = "red") +
    ggplot2::labs(x = "DA score (negative = depleted in Depletion)",
                  y = "Cells",
                  title = "DA-seq: score distribution with permutation thresholds") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  ggplot2::ggsave(plot = p,
                  filename = file.path(rd, "daseq_score_histogram.png"),
                  width = 6, height = 3, dpi = 300, bg = "white")
  rm(p)
}

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat("DA-seq cell-level DA detection on simulated data\n")
cat("Results saved to:", rd, "\n")
cat(rep("=", 60), "\n")

sessionInfo()
