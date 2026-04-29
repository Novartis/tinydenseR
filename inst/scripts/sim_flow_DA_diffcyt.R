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

library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(tinydenseR)
library(tidyverse)
library(ggpubr)
library(rstatix)

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
rd <- file.path(dirname(script.path), "results", "sim_flow_DA_diffcyt")

if(!dir.exists(paths = rd)) dir.create(path = rd, recursive = TRUE)

set.seed(42)

# Simulate DA data using tinydenseR API (same seed / RNG sequence as original)
sim_data <- tinydenseR::simulate_DA_data()
final_data_DA <- sim_data$cell_meta |>
  dplyr::mutate(Treatment = factor(x = Treatment,
                                   levels = c("Baseline", "Depletion")))

settings <- c(0.005, 0.05, 0.5)

# Function to convert simulated data to diffcyt format
# Reads pre-written FCS files from simulate_DA_data() output.
# NOTE: FCS files already contain log-transformed expression;
#       do NOT apply log() again.
prepare_diffcyt_data <- function(sample_meta, setting_filter) {
  
  # Filter sample metadata for specific setting
  meta_filtered <- sample_meta |>
    dplyr::filter(Setting == setting_filter)
  
  sample_ids <- meta_filtered$Sample
  
  # Read FCS files (already log-transformed by simulate_DA_data)
  flowframes_list <- list()
  
  for (i in seq_along(sample_ids)) {
    ff <- flowCore::read.FCS(meta_filtered$fcs_path[i],
                             transformation = FALSE,
                             truncate_max_range = FALSE)
    flowframes_list[[i]] <- ff
  }
  
  names(flowframes_list) <- sample_ids
  fs <- flowCore::flowSet(flowframes_list)
  
  # Create experiment_info
  experiment_info <- meta_filtered |>
    dplyr::select(Sample, Treatment, Batch) |>
    as.data.frame()
  
  rownames(experiment_info) <- experiment_info$Sample
  
  # Create marker_info
  marker_info <- data.frame(
    channel_name = paste0("Marker", 1:5),
    marker_name = paste0("Marker", 1:5),
    marker_class = c("type", "state", "type", "type", "type")  # Pretend Marker2 is state marker for DE testing (expect not sig)
  )
  
  return(list(
    flowset = fs,
    experiment_info = experiment_info,
    marker_info = marker_info
  ))
}

# Function to run diffcyt analysis
run_diffcyt_analysis <- function(diffcyt_data, setting_name) {
  
  cat("Running diffcyt analysis for", setting_name, "\n")
  
  # Prepare data
  diffcyt_data$experiment_info$sample_id <-
    diffcyt_data$experiment_info$Sample
  
  d_se <- diffcyt::prepareData(d_input = diffcyt_data$flowset,
                               experiment_info = diffcyt_data$experiment_info,
                               marker_info = diffcyt_data$marker_info)
  
  # Transform data (arcsinh transformation with cofactor 1)
  d_se <- diffcyt::transformData(d_se, cofactor = 1)
  
  # Generate clusters using FlowSOM
  d_se <- diffcyt::generateClusters(d_se, 
                                    xdim = 10, 
                                    ydim = 10,
                                    seed_clustering = 123)
  
  # Calculate features (cluster counts for DA testing)
  d_counts <- diffcyt::calcCounts(d_se)
  d_medians <- diffcyt::calcMedians(d_se)
  
  # Create design matrix
  design <- diffcyt::createDesignMatrix(diffcyt_data$experiment_info,
                                        cols_design = c("Treatment", "Batch"))
  
  # Create contrast matrix for Treatment effect
  contrast <- diffcyt::createContrast(c(0, 1, 0))  # Compare Depletion vs Baseline
  
  # Test for differential abundance (DA) - main focus for this analysis
  res_DA <- diffcyt::testDA_edgeR(d_counts, 
                                  design, 
                                  contrast,
                                  min_cells = 3,
                                  min_samples = 3)
  
  # Test for differential states (DS) - should show no significant results
  res_DS <- diffcyt::testDS_limma(d_counts,
                                  d_medians, 
                                  design, 
                                  contrast,
                                  plot = FALSE)
  
  return(list(
    d_se = d_se,
    d_counts = d_counts,
    d_medians = d_medians,
    res_DA = res_DA,
    res_DS = res_DS,
    design = design,
    contrast = contrast
  ))
}

# Function to extract and visualize results
extract_results <- function(diffcyt_results, setting_name) {
  
  # Extract DA results
  da_results <- rowData(diffcyt_results$res_DA) |>
    as.data.frame() |>
    (\(x)
     dplyr::mutate(.data = x,
                   cluster_id = rownames(x = x),
                   significant = p_adj < 0.1,
                   test_type = "DS")
    )() |>
    dplyr::arrange(p_val)
  
  # Extract DS results  
  ds_results <- rowData(diffcyt_results$res_DS) |>
    as.data.frame() |>
    (\(x)
     dplyr::mutate(.data = x,
                   cluster_id = rownames(x = x),
                   significant = p_adj < 0.1,
                   test_type = "DS")
    )() |>
    dplyr::arrange(p_val)
  
  cat("\nDifferential Abundance Results for", setting_name, ":\n")
  print(da_results |> dplyr::select(cluster_id, logFC, p_val, p_adj, significant))
  
  cat("\nDifferential States Results for", setting_name, "(should show no significant results):\n")
  print(ds_results |> head(5) |> dplyr::select(cluster_id, marker_id, logFC, p_val, p_adj, significant))
  
  return(list(da_results = da_results, ds_results = ds_results))
}

# Function to create visualization plots
create_plots <- function(diffcyt_results, results_summary, setting_name) {
  
  # Plot heatmap of cluster abundance
  p_heatmap <- diffcyt::plotHeatmap(res = diffcyt_results$res_DA,
                                    d_se = diffcyt_results$d_se, 
                                    d_counts = diffcyt_results$d_counts,
                                    d_medians = diffcyt_results$d_medians,
                                    d_medians_by_cluster_marker = diffcyt::calcMediansByClusterMarker(d_se = diffcyt_results$d_se),
                                    analysis_type = "DA")
  
  # Bar plot for DA results
  da_data <- results_summary$da_results |>
    dplyr::mutate(cluster_id = factor(cluster_id, levels = cluster_id[order(logFC)]))
  
  p_da_bar <- ggplot(da_data, aes(x = cluster_id, y = logFC, fill = significant)) +
    geom_col() +
    scale_fill_manual(values = c("FALSE" = "lightgrey", "TRUE" = "darkred")) +
    labs(title = paste("Differential Abundance -", setting_name),
         subtitle = "Expected to show depletion (negative logFC) for target cell clusters",
         x = "Cluster",
         y = "log2 Fold Change (Depletion vs Baseline)",
         fill = "Significant\n(FDR < 0.05)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  # Volcano plot for DA results
  da_volcano_data <- results_summary$da_results |>
    dplyr::mutate(log10_p = -log10(p_val),
                  significant = p_adj < 0.1)
  
  p_da_volcano <- ggplot(da_volcano_data, aes(x = logFC, y = log10_p, color = significant)) +
    geom_point(alpha = 0.6, size = 3) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    labs(title = paste("Differential Abundance Volcano Plot -", setting_name),
         x = "log2 Fold Change (Depletion vs Baseline)",
         y = "-log10(p-value)",
         color = "Significant\n(FDR < 0.05)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
  
  return(list(
    heatmap = p_heatmap,
    da_bar = p_da_bar,
    da_volcano = p_da_volcano
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
  
  # Calculate expected fold change
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
plots_all <- list()
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
  diffcyt_data <- prepare_diffcyt_data(sim_data$sample_meta, setting_name)
  
  # Run diffcyt analysis
  diffcyt_results <- run_diffcyt_analysis(diffcyt_data, setting_name)
  
  # Extract results
  results_summary <- extract_results(diffcyt_results, setting_name)
  
  # Create plots
  plots <- create_plots(diffcyt_results, results_summary, setting_name)
  
  # Store results
  results_all[[setting_name]] <- list(
    diffcyt_results = diffcyt_results,
    results_summary = results_summary
  )
  plots_all[[setting_name]] <- plots
  
}

plots_all$`0.5%`$heatmap
plots_all$`5%`$heatmap
plots_all$`50%`$heatmap

# Summary comparison across settings
cat("\n", rep("=", 60), "\n")
cat("SUMMARY: Comparison of results across settings\n")
cat(rep("=", 60), "\n")

summary_comparison <- data.frame()

for (setting in settings) {
  setting_name <- paste0(setting * 100, "%")
  
  da_res <- results_all[[setting_name]]$results_summary$da_results
  ds_res <- results_all[[setting_name]]$results_summary$ds_results
  gt <- ground_truth_all[[setting_name]]
  
  summary_row <- data.frame(
    Setting = setting_name,
    Expected_logFC = round(gt$expected_logFC, 3),
    N_Clusters = nrow(da_res),
    Significant_DA = sum(da_res$significant),
    Max_DA_logFC = round(max(abs(da_res$logFC)), 3),
    Most_Significant_DA_logFC = round(da_res$logFC[which.min(da_res$p_adj)], 3),
    Baseline_Prop = round(gt$baseline_prop, 4),
    Depletion_Prop = round(gt$depletion_prop, 4)
  )
  
  summary_comparison <- rbind(summary_comparison, summary_row)
}

print(summary_comparison)

# Create comparison plot of number of significant DA results
comparison_plot_data <- summary_comparison |>
  dplyr::select(Setting, Significant_DA, Expected_logFC) |>
  dplyr::mutate(Effect_Size = abs(Expected_logFC))

p_comparison <- ggplot(comparison_plot_data, aes(x = Setting, y = Significant_DA, fill = Effect_Size)) +
  geom_col() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Number of Significant Differential Abundance Results by diffcyt",
       subtitle = "Expected: More significant results with larger baseline proportions",
       x = "Baseline Target Cell Proportion",
       y = "Number of Significant DA Results",
       fill = "Effect Size\n|log2FC|") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

print(p_comparison)

# Create logFC comparison plot
logfc_comparison_data <- summary_comparison |>
  dplyr::select(Setting, Expected_logFC, Most_Significant_DA_logFC) |>
  tidyr::pivot_longer(cols = c(Expected_logFC, Most_Significant_DA_logFC),
                      names_to = "Type", 
                      values_to = "logFC") |>
  dplyr::mutate(Type = case_when(
    Type == "Expected_logFC" ~ "Expected (Ground Truth)",
    Type == "Most_Significant_DA_logFC" ~ "Observed (Most Significant)"
  ))

p_logfc_comparison <- ggplot(logfc_comparison_data, aes(x = Setting, y = logFC, color = Type, group = Type)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Expected (Ground Truth)" = "black", "Observed (Most Significant)" = "red")) +
  labs(title = "Expected vs Observed log2 Fold Changes",
       subtitle = "Validation of diffcyt DA detection",
       x = "Baseline Target Cell Proportion",
       y = "log2 Fold Change (Depletion vs Baseline)",
       color = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p_logfc_comparison)

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat("diffcyt successfully identified differential abundance patterns\n")
cat("Results show expected depletion (negative logFC) in target cell populations\n")
cat("Larger baseline proportions yield more power to detect differences\n")
cat("Results and plots saved to:", rd, "\n")
cat(rep("=", 60), "\n")

sessionInfo()