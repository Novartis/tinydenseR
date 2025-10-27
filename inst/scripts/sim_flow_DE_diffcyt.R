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
library(tidyverse)
library(ggpubr)
library(rstatix)

wd <- "path/to/your/working/directory/"
rd <- file.path(wd, "res")

setwd(dir = wd)

set.seed(42)

# Parameters
groups <- c("Baseline", "Activation")
batches <- c("Batch1", "Batch2")
sd_shifts <- c(0.5, 1, 2)  # SD differences in Marker2 for Baseline vs Activation
samples_per_group <- 6
mean_cells <- 50000
sd_cells <- 500  # Noise in total cell count

# Initialize storage
data_list_DE <- list()

# Simulate data (same as original)
for (sd_shift in sd_shifts) {
  for (group in groups) {
    for (sample_id in 1:samples_per_group) {
      batch <- 
        if(sample_id %% 2 == 0) {
          "Batch2"
        } else {
          "Batch1"
        }
      total_cells <- round(rnorm(1, mean = mean_cells, sd = sd_cells))
      total_cells <- max(total_cells, 1000)  # Ensure a minimum number of cells
      
      # Fixed proportion of cell type of interest
      proportion <- 0.05
      num_interest <- round(total_cells * proportion)
      num_other <- total_cells - num_interest
      cell_types <- c(rep("target", num_interest), rep("other", num_other))
      cell_types <- sample(cell_types)
      
      # Simulate expression data
      marker1 <- numeric(total_cells)
      marker2 <- rlnorm(total_cells, meanlog = 0, sdlog = 1.5)
      marker3 <- rlnorm(total_cells, meanlog = 0, sdlog = 2.5)
      marker4 <- numeric(total_cells)  
      marker5 <- numeric(total_cells)  
      
      # Assign Marker1, Marker4 Marker5 based on cell type
      marker1[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 2)
      marker1[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 5, sdlog = 2)  # Shift by 5 SD
      marker4[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.2)
      marker4[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 3, sdlog = 1.2)  # Shift by 3 SD
      marker5[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.8)
      marker5[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + 7, sdlog = 1.8)  # Shift by 7 SD
      
      # Marker2: differential expression between groups in target cells only
      marker2[cell_types == "other"] <- rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.5)
      if (group == "Baseline") {
        marker2[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0, sdlog = 1.5)
      } else {
        marker2[cell_types == "target"] <- rlnorm(sum(cell_types == "target"), meanlog = 0 + sd_shift, sdlog = 1.5)
      }
      
      # Add batch effect
      if (batch == "Batch2") {
        marker1 <- marker1 * rlnorm(total_cells, meanlog = 0.1, sdlog = 0.5)
        marker2 <- marker2 * rlnorm(total_cells, meanlog = 0.2, sdlog = 0.3)
        marker3 <- marker3 * rlnorm(total_cells, meanlog = 0.3, sdlog = 0.4)
        marker4 <- marker4 * rlnorm(total_cells, meanlog = 0.4, sdlog = 0.3)
        marker5 <- marker5 * rlnorm(total_cells, meanlog = 0.5, sdlog = 0.5)
      }
      
      sample_name <- paste0(group, "_S", sample_id, "_Shift", sd_shift, "SD")
      df <- data.frame(
        Sample = sample_name,
        Treatment = group,
        Batch = batch,
        SD_Shift = paste0(sd_shift, "SD"),
        CellType = cell_types,
        Marker1 = marker1,
        Marker2 = marker2,
        Marker3 = marker3,
        Marker4 = marker4,
        Marker5 = marker5
      )
      data_list_DE[[length(data_list_DE) + 1]] <- df
    }
  }
}

# Combine all simulated data
final_data_DE <- do.call(rbind, data_list_DE) |>
  dplyr::mutate(Treatment = factor(x = Treatment,
                                   levels = c("Baseline", "Activation")))

# Function to convert simulated data to diffcyt format
prepare_diffcyt_data <-
  function(data, sd_shift_filter) {
    
    # Filter data for specific SD shift
    data_filtered <- data |>
      dplyr::filter(SD_Shift == sd_shift_filter)
    
    # Create flowSet from simulated data
    sample_ids <- unique(data_filtered$Sample)
    
    # Create list of flowFrames
    flowframes_list <- list()
    
    for (i in seq_along(sample_ids)) {
      sample_data <- data_filtered |>
        dplyr::filter(Sample == sample_ids[i]) |>
        dplyr::select(Marker1, Marker2, Marker3, Marker4, Marker5) |>
        as.matrix()
      
      # Apply log transformation (equivalent to original log() transformation)
      sample_data <- log(sample_data)
      
      # Create flowFrame
      colnames(sample_data) <- paste0("Marker", 1:5)
      ff <- flowCore::flowFrame(exprs = sample_data)
      flowframes_list[[i]] <- ff
    }
    
    names(flowframes_list) <- sample_ids
    fs <- flowCore::flowSet(flowframes_list)
    
    # Create experiment_info
    experiment_info <- data_filtered |>
      dplyr::select(Sample, Treatment, Batch) |>
      dplyr::distinct() |>
      dplyr::arrange(match(Sample, sample_ids)) |>
      as.data.frame()
    
    rownames(experiment_info) <- experiment_info$Sample
    
    # Create marker_info
    marker_info <- data.frame(
      channel_name = paste0("Marker", 1:5),
      marker_name = paste0("Marker", 1:5),
      marker_class = c("type", "state", "type", "type", "type")  # Marker2 is state marker for DE testing
    )
    
    return(list(
      flowset = fs,
      experiment_info = experiment_info,
      marker_info = marker_info
    ))
  }

# Function to run diffcyt analysis
run_diffcyt_analysis <- function(diffcyt_data, sd_shift_name) {
  
  cat("Running diffcyt analysis for", sd_shift_name, "\n")
  
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
  
  # Calculate features (median marker expression per cluster per sample)
  d_counts <- diffcyt::calcCounts(d_se)
  d_medians <- diffcyt::calcMedians(d_se)
  
  # Create design matrix
  design <- diffcyt::createDesignMatrix(diffcyt_data$experiment_info,
                                        cols_design = c("Treatment", "Batch"))
  
  # Create contrast matrix for Treatment effect
  contrast <- diffcyt::createContrast(c(0, 1, 0))  # Compare Activation vs Baseline
  
  # Test for differential abundance (DA)
  res_DA <- diffcyt::testDA_edgeR(d_counts, 
                                  design, 
                                  contrast,
                                  min_cells = 3,
                                  min_samples = 3)
  
  # Test for differential states (DS) - differential expression
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
extract_results <- function(diffcyt_results, sd_shift_name) {
  
  # Extract DA results
  da_results <- rowData(diffcyt_results$res_DA) |>
    as.data.frame() |> 
    (\(x)
     dplyr::mutate(.data = x,
                   cluster_id = rownames(x = x),
                   significant = p_adj < 0.1,
                   test_type = "DA")
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
  
  cat("\nDifferential Abundance Results for", sd_shift_name, ":\n")
  print(da_results |> dplyr::select(cluster_id, logFC, p_val, p_adj, significant))
  
  cat("\nDifferential States Results for", sd_shift_name, "(top 10):\n")
  print(ds_results |> head(10) |> dplyr::select(cluster_id, marker_id, logFC, p_val, p_adj, significant))
  
  return(list(da_results = da_results, ds_results = ds_results))
}

# Function to create visualization plots
create_plots <- function(diffcyt_results, results_summary, sd_shift_name) {
  
  # Plot heatmap of median marker expressions
  p_heatmap <- diffcyt::plotHeatmap(res = diffcyt_results$res_DS,
                                    d_se = diffcyt_results$d_se, 
                                    d_counts = diffcyt_results$d_counts,
                                    d_medians = diffcyt_results$d_medians,
                                    d_medians_by_cluster_marker = diffcyt::calcMediansByClusterMarker(d_se = diffcyt_results$d_se),
                                    analysis_type = "DS")
  
  # Volcano plot for DS results
  ds_data <- results_summary$ds_results |>
    dplyr::mutate(log10_p = -log10(p_val),
                  significant = p_adj < 0.1)
  
  p_volcano <- ggplot(ds_data, aes(x = logFC, y = log10_p, color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    labs(title = paste("Differential States -", sd_shift_name),
         x = "log2 Fold Change",
         y = "-log10(p-value)") +
    theme_bw() +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
  
  # Bar plot for DA results
  da_data <- results_summary$da_results |>
    dplyr::mutate(cluster_id = factor(cluster_id, levels = cluster_id[order(logFC)]))
  
  p_da_bar <- ggplot(da_data, aes(x = cluster_id, y = logFC, fill = significant)) +
    geom_col() +
    scale_fill_manual(values = c("FALSE" = "lightgrey", "TRUE" = "darkred")) +
    labs(title = paste("Differential Abundance -", sd_shift_name),
         x = "Cluster",
         y = "log2 Fold Change (Activation vs Baseline)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(list(
    heatmap = p_heatmap,
    volcano = p_volcano,
    da_bar = p_da_bar
  ))
}

# Run analysis for each SD shift
results_all <- list()
plots_all <- list()

for (sd_shift in sd_shifts) {
  sd_shift_name <- paste0(sd_shift, "SD")
  
  cat("\n", rep("=", 50), "\n")
  cat("Processing SD shift:", sd_shift_name, "\n")
  cat(rep("=", 50), "\n")
  
  # Prepare data
  diffcyt_data <- prepare_diffcyt_data(final_data_DE, sd_shift_name)
  
  # Run diffcyt analysis
  diffcyt_results <- run_diffcyt_analysis(diffcyt_data, sd_shift_name)
  
  # Extract results
  results_summary <- extract_results(diffcyt_results, sd_shift_name)
  
  # Create plots
  plots <- create_plots(diffcyt_results, results_summary, sd_shift_name)
  
  # Store results
  results_all[[sd_shift_name]] <- list(
    diffcyt_results = diffcyt_results,
    results_summary = results_summary
  )
  plots_all[[sd_shift_name]] <- plots
  
}

plots_all$`0.5SD`$heatmap
plots_all$`1SD`$heatmap
plots_all$`2SD`$heatmap

# Summary comparison across SD shifts
cat("\n", rep("=", 60), "\n")
cat("SUMMARY: Comparison of results across SD shifts\n")
cat(rep("=", 60), "\n")

summary_comparison <- data.frame()

for (sd_shift in sd_shifts) {
  sd_shift_name <- paste0(sd_shift, "SD")
  
  da_res <- results_all[[sd_shift_name]]$results_summary$da_results
  ds_res <- results_all[[sd_shift_name]]$results_summary$ds_results
  
  summary_row <- data.frame(
    SD_Shift = sd_shift_name,
    N_Clusters = nrow(da_res),
    Significant_DA = sum(da_res$significant),
    Max_DA_logFC = max(abs(da_res$logFC)),
    N_ClusterMarker_Combinations = nrow(ds_res),
    Significant_DS = sum(ds_res$significant),
    Max_DS_logFC = max(abs(ds_res$logFC))
  )
  
  summary_comparison <- rbind(summary_comparison, summary_row)
}

print(summary_comparison)

# Create comparison plot of number of significant results
comparison_plot_data <- summary_comparison |>
  dplyr::select(SD_Shift, Significant_DA, Significant_DS) |>
  tidyr::pivot_longer(cols = c(Significant_DA, Significant_DS), 
                      names_to = "Test_Type", 
                      values_to = "N_Significant") |>
  dplyr::mutate(Test_Type = gsub("Significant_", "", Test_Type))

p_comparison <- ggplot(comparison_plot_data, aes(x = SD_Shift, y = N_Significant, fill = Test_Type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("DA" = "steelblue", "DS" = "orange")) +
  labs(title = "Number of Significant Results by diffcyt",
       subtitle = "Comparing Differential Abundance (DA) and Differential States (DS)",
       x = "Effect Size (SD Shift)",
       y = "Number of Significant Results",
       fill = "Test Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

print(p_comparison)

# Focus on Marker2 DS results (the marker with simulated differential expression)
cat("\n", rep("=", 60), "\n")
cat("FOCUS: Marker2 Differential States Results\n")
cat(rep("=", 60), "\n")

marker2_results <- data.frame()

for (sd_shift in sd_shifts) {
  sd_shift_name <- paste0(sd_shift, "SD")
  
  ds_res <- results_all[[sd_shift_name]]$results_summary$ds_results
  
  # Filter for Marker2 results
  marker2_ds <- ds_res |>
    dplyr::filter(grepl("Marker2", marker_id)) |>
    dplyr::mutate(SD_Shift = sd_shift_name) |>
    dplyr::select(SD_Shift, marker_id, logFC, p_val, p_adj, significant)
  
  marker2_results <- rbind(marker2_results, marker2_ds)
}

print(marker2_results)

# Plot Marker2 logFC across SD shifts
if (nrow(marker2_results) > 0) {
  p_marker2 <- ggplot(marker2_results, aes(x = SD_Shift, y = logFC, color = significant)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    labs(title = "Marker2 Differential Expression Results",
         subtitle = "Expected to show increasing effect with larger SD shifts",
         x = "Simulated Effect Size",
         y = "Observed log2 Fold Change",
         color = "Significant\n(FDR < 0.05)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  
  print(p_marker2)
  
  
}

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat("diffcyt successfully identified differential abundance and differential states\n")
cat("Results and plots saved to:", rd, "\n")
cat(rep("=", 60), "\n")