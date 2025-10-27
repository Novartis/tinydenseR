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

library(tinydenseR)
library(Seurat)
library(SeuratObject)
library(BPCells)
library(Azimuth)
library(miloR)
library(SingleCellExperiment)
library(BiocSingular)
library(DelayedArray)
library(DelayedMatrixStats)
library(HDF5Array)
library(biomaRt)
library(tidyverse)
library(profvis)
library(ggpubr)
options(future.globals.maxSize = 1e10)

wd <- 
  "path/to/your/working/directory" # change to your working directory
rd <-
  file.path(wd,
            "res")

setwd(dir = wd)

# Load the peak_mem_profvis (see inst/reproducible_analysis_script/) function
# for accurate memory and timing measurements
source(file = file.path("path/to/your/file/location",
                        "peak_mem_profvis.R"))

# Follow instructions to create Seurat object from h5ad files here
# https://github.com/satijalab/seurat/blob/e35abd442520808a20025e589f861620ddc315af/vignettes/seurat5_bpcells_interaction_vignette.Rmd
# and just to line 59 here: 
# https://github.com/satijalab/seurat/blob/30f82df52159ac5f0feb80b149698abbd876b779/vignettes/COVID_SCTMapping.Rmd
# and save the object to the disk
processed_obj_path <- 
  file.path("path/to/your/file/location",
            "object.Rds")

saveRDS(object = object,
        file = processed_obj_path)

rm(object)

# Set up for reproducibility
set.seed(42)

# Parameters for benchmarking
cell_counts <- 2^(14:21)  # 8 different cell count levels: 16K, 32K, 65K, 131K, 262K, 524K, 1M, all
n_runs <- 3  # 3 runs per method per cell count (63 total analyses)

# Initialize results storage
benchmark_results <- data.frame()

# Function to load and prepare data
load_data <- function() {
  
  cat("Loading and preparing COVID dataset...\n")
  
  # Load reference data
  reference <- 
    file.path("path/to/your/file/location",
              "pbmc_multimodal_2023.rds") |>
    readRDS()
  
  # Check if processed object already exists
  processed_obj_path <- 
    file.path("path/to/your/file/location",
              "object.Rds") |>
    readRDS()
  
  cat("Loading pre-processed merged object...\n")
  merged.object <- readRDS(file = processed_obj_path)
  
  cat("Data loaded successfully. Total cells:", ncol(merged.object), "\n")
  cat("Total genes:", nrow(merged.object), "\n")
  cat("Publications:", unique(merged.object$publication), "\n")
  cat("Diseases:", unique(merged.object$disease), "\n")
  
  return(list(reference = reference, merged = merged.object))
}

# Function to subset data to specific cell count
subset_data <- function(merged_object, n_cells) {
  
  # Get total available cells
  total_cells <- ncol(merged_object)
  
  cat(paste("Targeting", format(n_cells, big.mark = ","), "cells from", format(total_cells, big.mark = ","), "total cells\n"))
  
  # Get metadata and create sample summary
  metadata <- merged_object@meta.data
  
  if("donor_id_disease" %in% colnames(metadata)) {
    cat("donor_id_disease column already exists in metadata.\n")
  } else {
    cat("Adding donor_id_disease column to metadata.\n")
    # Add donor_id_disease column (as in original script)
    merged_object$donor_id_disease <- paste0(merged_object$donor_id, "_", merged_object$disease) |>
      make.names()
    # Update metadata after adding the column
    metadata <- merged_object@meta.data
  }
  
  # Use only samples with more than a thousand cells to reduce noise, see:
  # table(merged_object$donor_id_disease) |>
  #  sort() |>
  #  plot(ylab = "count"); abline(a = 1000, b = 0)
  valid_samples <-
    table(merged_object$donor_id_disease) |> 
    (\(x)
     x[x > 1000]
    )() |>
    names()
  
  metadata <- metadata[metadata$donor_id_disease %in% valid_samples,]
  
  # Create sample summary: count cells per donor_id_disease and get publication info
  sample_summary <- metadata %>%
    dplyr::group_by(donor_id_disease, publication, disease) %>%
    dplyr::summarise(
      n_cells_in_sample = n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(n_cells_in_sample))  # Start with larger samples for efficiency
  
  cat("Sample summary:\n")
  cat(paste("Total samples available:", nrow(sample_summary), "\n"))
  cat("Cells per publication:\n")
  print(sample_summary %>% 
          dplyr::group_by(publication) %>% 
          dplyr::summarise(samples = n(), total_cells = sum(n_cells_in_sample), .groups = "drop"))
  
  # Sample selection strategy: maintain publication balance while approaching target
  selected_samples <- c()
  current_cell_count <- 0
  
  # Calculate target samples per publication (proportional to their representation)
  pub_proportions <- sample_summary %>%
    dplyr::group_by(publication) %>%
    dplyr::summarise(
      pub_cells = sum(n_cells_in_sample),
      pub_samples = n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      proportion = pub_samples / sum(pub_samples),
      target_samples = pmax(1, round(proportion * (n_cells / mean(sample_summary$n_cells_in_sample))))
    )
  
  cat("Target samples per publication:\n")
  print(pub_proportions)
  
  # Sample from each publication proportionally
  set.seed(42)  # Ensure reproducible sampling
  
  for (pub in pub_proportions$publication) {
    pub_samples <- sample_summary %>% dplyr::filter(publication == pub)
    target_for_pub <- pub_proportions$target_samples[pub_proportions$publication == pub]
    
    # Don't exceed available samples for this publication
    n_to_sample <- min(target_for_pub, nrow(pub_samples))
    
    if (n_to_sample > 0) {
      # Randomly sample from this publication
      selected_pub_samples <- sample(pub_samples$donor_id_disease, n_to_sample)
      selected_samples <- c(selected_samples, selected_pub_samples)
      
      # Update cell count
      pub_cell_count <- sum(pub_samples$n_cells_in_sample[pub_samples$donor_id_disease %in% selected_pub_samples])
      current_cell_count <- current_cell_count + pub_cell_count
      
      cat(paste("Selected", n_to_sample, "samples from", pub, "contributing", 
                format(pub_cell_count, big.mark = ","), "cells\n"))
    }
  }
  
  # If we're still significantly under target, add more samples from any publication
  if ((current_cell_count < n_cells * 0.8) & # If we're more than 20% under target
      (n_cells < 2^21)) {  # Avoid overshooting
    remaining_samples <- sample_summary %>%
      dplyr::filter(!donor_id_disease %in% selected_samples) %>%
      dplyr::arrange(desc(n_cells_in_sample))  # Prefer larger samples
    
    for (i in 1:nrow(remaining_samples)) {
      if (current_cell_count >= n_cells * 0.8) break  # Stop when we're close enough
      
      additional_sample <- remaining_samples$donor_id_disease[i]
      additional_cells <- remaining_samples$n_cells_in_sample[i]
      
      selected_samples <- c(selected_samples, additional_sample)
      current_cell_count <- current_cell_count + additional_cells
      
      cat(paste("Added sample", additional_sample, "with", 
                format(additional_cells, big.mark = ","), "cells\n"))
    }
  }
  
  # Get all cells from selected samples
  selected_cells <- rownames(metadata)[metadata$donor_id_disease %in% selected_samples]
  
  # Extract components efficiently (much faster than Seurat subsetting with BPCells)
  cat("Extracting count matrix and metadata...\n")
  
  # for each layer, get counts from all cells from selected samples
  # get intersect of genes
  .genes.names <-
    names(x = merged_object@assays$RNA@layers) |>
    grep(pattern = "counts",
         value = TRUE,
         fixed = TRUE) |>
    lapply(FUN = function(dataset){
      
      SeuratObject::LayerData(object = merged_object, 
                              assay = "RNA", 
                              layer = dataset) |>
        rownames()
    }) |>
    Reduce(f = base::intersect)
  
  subset_counts <-
    names(x = merged_object@assays$RNA@layers) |>
    grep(pattern = "counts",
         value = TRUE,
         fixed = TRUE) |>
    lapply(FUN = function(data.layer){
      
      original_counts_layer <-
        SeuratObject::LayerData(object = merged_object, 
                                assay = "RNA", 
                                layer = data.layer)
      
      return(original_counts_layer[.genes.names,
                                   colnames(x = original_counts_layer) %in%
                                     selected_cells])
    }) |>
    do.call(what = cbind) |> 
    (\(x)
     x[,selected_cells]
    )()
  
  # Preserve BPCells format by subsetting the BPCells matrix directly
  # Check if it's a BPCells matrix and preserve the format
  if (inherits(subset_counts, c("IterableMatrix", "RealizationSink"))) {
    cat("Preserving BPCells matrix format...\n") 
  } else {
    stop("Matrix is not BPCells format, using standard subsetting...\n")
  }
  
  subset_metadata <- merged_object@meta.data[selected_cells, ]
  
  # Create new Seurat object (much faster than subsetting)
  cat("Creating new Seurat object...\n")
  subset_obj <- Seurat::CreateSeuratObject(
    counts = subset_counts, 
    meta.data = subset_metadata,
    project = "benchmark_subset"
  )
  
  # Verify BPCells format is preserved
  if (inherits(subset_obj@assays$RNA@layers$counts, c("IterableMatrix", "RealizationSink"))) {
    cat("✓ BPCells format preserved in new Seurat object\n")
  } else {
    stop("⚠ BPCells format was not preserved\n")
  }
  
  actual_cells <- ncol(subset_obj)
  cat(paste("Successfully selected", length(selected_samples), "samples with", 
            format(actual_cells, big.mark = ","), "cells\n"))
  cat(paste("Target:", format(n_cells, big.mark = ","), "| Actual:", 
            format(actual_cells, big.mark = ","), 
            "| Difference:", format(actual_cells - n_cells, big.mark = ","), "\n"))
  
  cat("Final sample composition:\n")
  print(table(subset_obj$publication))
  cat("Final disease composition:\n") 
  print(table(subset_obj$disease))
  
  return(subset_obj)
}

# Function to benchmark Seurat workflow
benchmark_seurat <- function(query_obj, reference_obj, run_id, n_cells) {
  
  cat(paste("Running Seurat benchmark - Run", run_id, "with", format(n_cells, big.mark = ","), "cells\n"))
  
  # Use profvis to capture real-time memory usage (following original approach)
  prof_result <- 
    withCallingHandlers(
      profvis::profvis(expr = {
        # Normalize data
        query_normalized <- Seurat::NormalizeData(object = query_obj, verbose = FALSE)
        
        # Find transfer anchors with parameters optimized for large datasets
        anchors <- Seurat::FindTransferAnchors(
          reference = reference_obj,
          query = query_normalized,
          reference.reduction = "spca",
          normalization.method = "SCT"
        )
        
        # Map query to reference
        query_mapped <- Seurat::MapQuery(
          anchorset = anchors,
          query = query_normalized,
          reference = reference_obj,
          refdata = list(
            celltype.l1 = "celltype.l1",
            celltype.l2 = "celltype.l2"
          ),
          reduction.model = "wnn.umap"
        )
        
        return(query_mapped)
      }),
      message = function(m) {
        msg <- conditionMessage(m)
        if (startsWith(msg, "profvis: code exited with error:")) stop(msg)
      })
  
  # Extract timing and memory using peak_mem_profvis (more accurate than manual calculation)
  profvis_results <- peak_mem_profvis(prof_result, input = "auto", output = "MiB", verbose = FALSE)
  
  # Extract timing and memory information
  result <- data.frame(
    method = "Seurat",
    run_id = run_id,
    n_cells = n_cells,
    time_seconds = profvis_results$timing$elapsed_s,
    memory_mib = profvis_results$peak,
    peak_memory_bytes = profvis_results$peak_bytes,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Function to benchmark tinydenseR workflow
benchmark_tinydenseR <- function(query_obj, reference_obj, run_id, n_cells) {
  
  cat(paste("Running tinydenseR benchmark - Run", run_id, "with", format(n_cells, big.mark = ","), "cells\n"))
  
  # Prepare data for tinydenseR
  # Extract expression matrix and metadata
  expr_matrix <- SeuratObject::LayerData(query_obj, assay = "RNA", layer = "counts")
  metadata <- query_obj@meta.data
  
  # Create cells list for tinydenseR
  unique_samples <- unique(metadata$donor_id_disease)
  cells_list <- list()
  
  for (sample in unique_samples) {
    sample_cells <- rownames(metadata)[metadata$donor_id_disease == sample]
    # Save matrix to temporary file
    temp_file <- tempfile(fileext = ".RDS")
    sample_matrix <- as(object = expr_matrix[, sample_cells],
                        Class = "dgCMatrix")
    saveRDS(sample_matrix, temp_file, compress = FALSE)
    cells_list[[sample]] <- temp_file
    
  }
  
  cat(paste("Selected", length(cells_list), "samples for tinydenseR analysis\n"))
  
  if (length(cells_list) < 3) {
    stop("Not enough samples with sufficient cells for tinydenseR analysis")
  }
  
  # Create metadata for samples
  sample_metadata <- metadata[!duplicated(metadata$donor_id_disease), ]
  rownames(sample_metadata) <- sample_metadata$donor_id_disease
  sample_metadata <- sample_metadata[names(cells_list),]
  
  # Use profvis to capture real-time memory usage (following original approach)
  prof_result <- 
    withCallingHandlers(
      profvis::profvis(expr = {
        
        # Setup landmark object with optimized parameters for large datasets
        lm_obj <- tinydenseR::setup.lm.obj(
          .cells = cells_list,
          .meta = sample_metadata,
          .harmony.var = "publication",
          .assay.type = "RNA",
          .verbose = TRUE
        )
        
        # Get landmarks with parameters scaled for dataset size
        lm_obj <- tinydenseR::get.landmarks(
          .lm.obj = lm_obj, 
          .nHVG = 2000,
          .verbose = TRUE
        )
        
        # Get landmark graph with appropriate resolution
        lm_obj <- tinydenseR::get.graph(
          .lm.obj = lm_obj,
          .verbose = TRUE
        )
        
        # Map to reference
        lm_obj <- tinydenseR::get.map(
          .lm.obj = lm_obj,
          .ref.obj = reference_obj,
          .integrate.vars = "LibraryId",
          .celltype.col.name = "celltype.l2",
          .verbose = TRUE
        )
        return(lm_obj)
      }),
      message = function(m) {
        msg <- conditionMessage(m)
        if (startsWith(msg, "profvis: code exited with error:")) stop(msg)
      })
  
  # Clean up temp files
  sapply(cells_list, function(x) if(file.exists(x)) file.remove(x))
  
  # Extract timing and memory using peak_mem_profvis (more accurate than manual calculation)
  profvis_results <- peak_mem_profvis(prof_result, input = "auto", output = "MiB", verbose = FALSE)
  
  # Extract timing and memory information
  result <- data.frame(
    method = "tinydenseR",
    run_id = run_id,
    n_cells = n_cells,
    time_seconds = profvis_results$timing$elapsed_s,
    memory_mib = profvis_results$peak,
    peak_memory_bytes = profvis_results$peak_bytes,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Function to benchmark miloR with DelayedMatrix workflow
benchmark_miloR <- function(query_obj, run_id, n_cells) {
  
  cat(paste("Running miloR benchmark - Run", run_id, "with", format(n_cells, big.mark = ","), "cells\n"))
  
  # Convert Seurat to SingleCellExperiment with DelayedArray backend
  # Extract count matrix and convert to DelayedArray
  count_matrix <- SeuratObject::LayerData(query_obj, assay = "RNA", layer = "counts") |>
    as(Class = "dgCMatrix") |>
    DelayedArray::DelayedArray()
  
  # Use profvis to capture real-time memory usage (following original approach)
  prof_result <- 
    withCallingHandlers(
      profvis::profvis(expr = {
        
        # Create SingleCellExperiment
        sce <- SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = count_matrix),
          colData = query_obj@meta.data
        )
        
        libsizes <- Matrix::colSums(count_matrix)
        size.factors <- libsizes/mean(libsizes)
        SingleCellExperiment::logcounts(object = sce) <- log2(t(t(count_matrix)/size.factors) + 1) |>
          DelayedArray::DelayedArray()
        
        HVG <- 
          SingleCellExperiment::logcounts(object = sce) |>
          DelayedMatrixStats::rowVars() |>
          sort(decreasing = TRUE) |> 
          (\(x)
           names(x = x[1:2000])
          )() 
        
        SingleCellExperiment::reducedDim(x = sce, 
                                         type = "PCA") <-
          SingleCellExperiment::logcounts(object = sce)[HVG,] |>
          t() |>
          (\(x)
           BiocSingular::runPCA(x = x,
                                rank = 30,
                                center = DelayedMatrixStats::colMeans2(x),
                                scale = DelayedMatrixStats::colSds(x), 
                                BSPARAM = BiocSingular::IrlbaParam(deferred=TRUE))$x
          )()
        
        # Create Milo object
        milo_obj <- miloR::Milo(sce)
        
        # Build KNN graph
        milo_obj <- miloR::buildGraph(x = milo_obj, d = 30,
                                      transposed = TRUE)
        
        # Make neighborhoods
        milo_obj <- miloR::makeNhoods(milo_obj, refined = TRUE, )
        
        # Count cells in neighborhoods
        milo_obj <- miloR::countCells(milo_obj,
                                      meta.data = colData(milo_obj),
                                      samples = "donor_id_disease")
        
        return(milo_obj)
      }),
      message = function(m) {
        msg <- conditionMessage(m)
        if (startsWith(msg, "profvis: code exited with error:")) stop(msg)
      })
  
  # Extract timing and memory using peak_mem_profvis (more accurate than manual calculation)
  profvis_results <- peak_mem_profvis(prof_result, input = "auto", output = "MiB", verbose = FALSE)
  
  # Extract timing and memory information
  result <- data.frame(
    method = "miloR_DelayedMatrix",
    run_id = run_id,
    n_cells = n_cells,
    time_seconds = profvis_results$timing$elapsed_s,
    memory_mib = profvis_results$peak,
    peak_memory_bytes = profvis_results$peak_bytes,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Main benchmarking function
run_benchmarks <- function() {
  
  cat("Starting comprehensive benchmarking...\n")
  
  # Load data
  data_objects <- load_data()
  reference_obj <- data_objects$reference
  merged_obj <- data_objects$merged
  
  tdr.reference <-
    readRDS(file = "/da/onc/cBIT/milanpe1/dev/projects/cyto/manuscript/Seurat.data/symphony_ref.RDS")
  
  # Initialize results
  all_results <- data.frame()
  
  # Loop through each cell count
  for (n_cells in cell_counts) {
    
    cat("\n", rep("=", 60), "\n")
    cat("BENCHMARKING WITH", n_cells, "CELLS\n")
    cat(rep("=", 60), "\n")
    
    # Subset data for this cell count
    subset_obj <- subset_data(merged_obj, n_cells)
    actual_cells <- ncol(subset_obj)
    
    cat("Actual cells in subset:", actual_cells, "\n")
    
    # Run benchmarks for each method and each run
    for (run_id in 1:n_runs) {
      
      cat("\n--- Run", run_id, "of", n_runs, "---\n")
      
      # Seurat benchmark
      seurat_result <- benchmark_seurat(query_obj = subset_obj, 
                                        reference_obj = reference_obj, 
                                        run_id = run_id, 
                                        n_cells = actual_cells)
      all_results <- rbind(all_results, seurat_result)
      
      # Force garbage collection between runs
      gc()
      
      # tinydenseR benchmark  
      tinydense_result <- benchmark_tinydenseR(subset_obj, tdr.reference, run_id, actual_cells)
      all_results <- rbind(all_results, tinydense_result)
      
      # Force garbage collection between runs
      gc()
      
      # miloR benchmark
      if(actual_cells < 262144){
        milo_result <- benchmark_miloR(subset_obj, run_id, actual_cells)
        all_results <- rbind(all_results, milo_result)
        # Force garbage collection between runs
        gc()
      }
      
      # Save intermediate results
      saveRDS(all_results, file.path(rd, "benchmark_results_intermediate.rds"))
    }
  }
  
  return(all_results)
}

# Function to analyze and visualize results
analyze_results <- function(results) {
  
  cat("Analyzing benchmark results...\n")
  
  # Calculate summary statistics
  summary_stats <- results %>%
    group_by(method, n_cells) %>%
    summarise(
      mean_time = mean(time_seconds),
      sd_time = sd(time_seconds),
      mean_memory = mean(memory_mib),
      sd_memory = sd(memory_mib),
      n_runs = n(),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  # Create visualizations
  
  # 1. Execution time comparison
  p_time <- ggplot(results, aes(x = n_cells, y = time_seconds, color = method)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "Execution Time Comparison",
      subtitle = "Reference mapping and differential analysis across methods",
      x = "Number of Cells",
      y = "Execution Time (seconds)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # 2. Memory usage comparison
  p_memory <- ggplot(results, aes(x = n_cells, y = memory_mib, color = method)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "Memory Usage Comparison",
      subtitle = "Peak memory allocation during analysis",
      x = "Number of Cells",
      y = "Memory Usage (MiB)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # 3. Efficiency comparison (cells per second)
  results$efficiency <- results$n_cells / results$time_seconds
  
  p_efficiency <- ggplot(results, aes(x = n_cells, y = efficiency, color = method)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(
      title = "Processing Efficiency Comparison",
      subtitle = "Cells processed per second",
      x = "Number of Cells",
      y = "Efficiency (cells/second)",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # 4. Summary boxplots
  p_time_box <- ggplot(results, aes(x = method, y = time_seconds, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_log10() +
    labs(
      title = "Execution Time Distribution",
      x = "Method",
      y = "Time (seconds)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  p_memory_box <- ggplot(results, aes(x = method, y = memory_mib, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_log10() +
    labs(
      title = "Memory Usage Distribution",
      x = "Method",
      y = "Memory (MiB)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Save summary table
  write.csv(summary_stats, file.path(rd, "benchmark_summary_stats.csv"), row.names = FALSE)
  
  return(list(
    summary_stats = summary_stats,
    plots = list(
      time = p_time,
      memory = p_memory,
      efficiency = p_efficiency,
      time_box = p_time_box,
      memory_box = p_memory_box
    )
  ))
}

# Run the complete benchmarking pipeline
cat("Starting comprehensive benchmarking of tinydenseR, Seurat, and miloR...\n")
cat("Test Parameters:\n")
cat("- Cell counts: 100K, 250K, 500K, 1M, 1.5M cells\n")
cat("- Replicates: 3 runs per method per cell count (45 total analyses)\n")
cat("- Dataset: COVID PBMC data (ahern, jin, yoshida publications)\n")
cat("- Metrics: Execution time, memory usage, processing efficiency\n")
cat("This may take several hours to complete.\n\n")

# Run benchmarks
benchmark_results <- run_benchmarks()

# Save complete results
saveRDS(benchmark_results, file.path(rd, "benchmark_results_complete.rds"))
write.csv(benchmark_results, file.path(rd, "benchmark_results_complete.csv"), row.names = FALSE)

# Analyze results
analysis_results <- analyze_results(benchmark_results)

cat("\n", rep("=", 60), "\n")
cat("BENCHMARKING COMPLETE!\n")
cat(rep("=", 60), "\n")
cat("Results saved to:\n")
cat("- Raw data:", file.path(rd, "benchmark_results_complete.rds"), "\n")
cat("- CSV file:", file.path(rd, "benchmark_results_complete.csv"), "\n")
cat("- Summary stats:", file.path(rd, "benchmark_summary_stats.csv"), "\n")
cat("- Trend plots:", file.path(rd, "benchmark_comparison_trends.png"), "\n")
cat("- Distribution plots:", file.path(rd, "benchmark_comparison_distributions.png"), "\n")
cat(rep("=", 60), "\n")

# Print final summary
cat("\nFINAL SUMMARY:\n")
print(analysis_results$summary_stats)

cat("\nBenchmarking pipeline completed successfully!\n")

benchmark_results_complete <-
  readRDS(file = file.path(rd, "benchmark_results_complete.rds"))

benchmark_results_complete |>
  dplyr::mutate(method = ifelse(test = method == "Seurat",
                                yes = "Seurat_BPCells",
                                no = method)) |> 
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = n_cells,
                                          y = time_seconds/60,
                                          color = method,
                                          group = method)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(hjust = 0.5),
                    legend.title = ggplot2::element_blank()) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 3)))) +
     ggplot2::labs(
       title = "Processing Time",
       x = "number of cells",
       y = "time (minutes)"
     ) +
     ggplot2::geom_line(data = dplyr::group_by(.data = x, method, n_cells) |>
                          dplyr::summarize(time_seconds = median(time_seconds),
                                           .groups = "drop"),
                        linewidth = I(x = 0.25)) +
     ggplot2::scale_x_continuous(labels = c("0", "500K", "1M", "1.5M"),
                                 breaks = c(0, 5e5, 1e6, 1.5e6)) +
     ggplot2::geom_point(size = I(x = 1)) +
     ggplot2::scale_color_manual(
       values = grDevices::colorRampPalette(
         colors = unname(
           obj = tinydenseR::Color.Palette[1,c(1,7,2)])
       )(length(x = unique(x = x$method)))) +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 2,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 1.5,
                                                  units = "in"))
  )()

benchmark_results_complete |>
  dplyr::mutate(method = ifelse(test = method == "Seurat",
                                yes = "Seurat_BPCells",
                                no = method)) |> 
  (\(x)
   ggplot2::ggplot(data = x,
                   mapping = ggplot2::aes(x = n_cells,
                                          y = memory_mib * (2^20 / 10^9),
                                          color = method,
                                          group = method)) +
     ggplot2::theme_bw() +
     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(hjust = 0.5),
                    legend.title = ggplot2::element_blank()) +
     ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = I(x = 3)))) +
     ggplot2::labs(
       title = "Memory Usage",
       x = "number of cells",
       y = "memory (GB)"
     ) +
     ggplot2::geom_line(data = dplyr::group_by(.data = x, method, n_cells) |>
                          dplyr::summarize(memory_mib = median(memory_mib),
                                           .groups = "drop"),
                        linewidth = I(x = 0.25)) +
     ggplot2::scale_x_continuous(labels = c("0", "500K", "1M", "1.5M"),
                                 breaks = c(0, 5e5, 1e6, 1.5e6)) +
     ggplot2::geom_point(size = I(x = 1)) +
     ggplot2::scale_color_manual(
       values = grDevices::colorRampPalette(
         colors = unname(
           obj = tinydenseR::Color.Palette[1,c(1,7,2)])
       )(length(x = unique(x = x$method)))) +
     ggh4x::force_panelsizes(cols = ggplot2::unit(x = 2,
                                                  units = "in"),
                             rows = ggplot2::unit(x = 1.5,
                                                  units = "in"))
  )()
