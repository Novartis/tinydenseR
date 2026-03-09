#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

#' Simulate differential abundance (DA) flow cytometry data
#'
#' Generates simulated flow cytometry data with differential abundance of a
#' target cell type across treatment groups, batch effects, and multiple
#' abundance settings. Writes per-sample FCS files and returns sample- and
#' cell-level metadata.
#'
#' @param groups Character vector of treatment group names.
#' @param batches Character vector of batch names.
#' @param settings Numeric vector of proportions for the Baseline group.
#'   The non-Baseline group uses half the proportion.
#' @param samples_per_group Integer. Number of samples per group per setting.
#' @param mean_cells Integer. Mean number of cells per sample.
#' @param sd_cells Numeric. Standard deviation of cell count.
#' @param seed Integer. Random seed for reproducibility.
#' @param output_dir Character. Directory path where FCS files are written.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{sample_meta}{A \code{data.frame} with one row per sample
#'       and columns \code{Sample}, \code{Treatment}, \code{Batch},
#'       \code{Setting}, \code{fcs_path}.}
#'     \item{cell_meta}{A \code{data.frame} with one row per cell
#'       and columns \code{Sample}, \code{Treatment}, \code{Batch},
#'       \code{Setting}, \code{CellType}.}
#'   }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_DA_data(
#'   samples_per_group = 2,
#'   mean_cells = 1000,
#'   sd_cells = 100,
#'   output_dir = file.path(tempdir(), "sim_DA_fcs")
#' )
#' head(sim$sample_meta)
#' head(sim$cell_meta)
#' }
#'
#' @export
simulate_DA_data <- function(groups = c("Baseline", "Depletion"),
                             batches = c("Batch1", "Batch2"),
                             settings = c(0.005, 0.05, 0.5),
                             samples_per_group = 6,
                             mean_cells = 50000,
                             sd_cells = 500,
                             seed = 42,
                             output_dir = file.path(tempdir(), "sim_DA_fcs")) {

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required for simulate_DA_data(). ",
         "Install it with: BiocManager::install('flowCore')", call. = FALSE)
  }

  set.seed(seed)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize storage
  data_list <- list()

  # Simulate data
  for (setting in settings) {
    for (group in groups) {
      for (sample_id in 1:samples_per_group) {
        batch <-
          if (sample_id %% 2 == 0) {
            "Batch2"
          } else {
            "Batch1"
          }
        total_cells <- round(stats::rnorm(1, mean = mean_cells, sd = sd_cells))
        total_cells <- max(total_cells, 1000)  # Ensure a minimum number of cells

        # Two-fold difference in proportions
        proportion <- if (group == "Baseline") setting else setting / 2

        num_interest <- round(total_cells * proportion)
        num_other <- total_cells - num_interest
        cell_types <- c(rep("target", num_interest), rep("other", num_other))
        cell_types <- sample(cell_types)

        # Simulate expression data
        marker1 <- numeric(total_cells)
        marker2 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 1.5)
        marker3 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 2.5)
        marker4 <- numeric(total_cells)
        marker5 <- numeric(total_cells)

        # Assign Marker1, Marker4 Marker5 based on cell type
        marker1[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 2)
        marker1[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 5, sdlog = 2)  # Shift by 5 SD
        marker4[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.2)
        marker4[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 3, sdlog = 1.2)  # Shift by 3 SD
        marker5[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.8)
        marker5[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 7, sdlog = 1.8)  # Shift by 7 SD

        # Add batch effect
        if (batch == "Batch2") {
          marker1 <- marker1 * stats::rlnorm(total_cells, meanlog = 0.1, sdlog = 0.5)
          marker2 <- marker2 * stats::rlnorm(total_cells, meanlog = 0.2, sdlog = 0.3)
          marker3 <- marker3 * stats::rlnorm(total_cells, meanlog = 0.3, sdlog = 0.4)
          marker4 <- marker4 * stats::rlnorm(total_cells, meanlog = 0.4, sdlog = 0.3)
          marker5 <- marker5 * stats::rlnorm(total_cells, meanlog = 0.5, sdlog = 0.5)
        }

        sample_name <- paste0(group, "_S", sample_id, "_Set", setting * 100)
        df <- data.frame(
          Sample = sample_name,
          Treatment = group,
          Batch = batch,
          Setting = paste0(setting * 100, "%"),
          CellType = cell_types,
          Marker1 = marker1,
          Marker2 = marker2,
          Marker3 = marker3,
          Marker4 = marker4,
          Marker5 = marker5
        )
        data_list[[length(data_list) + 1]] <- df
      }
    }
  }

  # Combine
  final_data <- do.call(rbind, data_list)

  # Write FCS files and build sample_meta
  sample_names <- unique(final_data$Sample)
  sample_meta_list <- list()

  for (sample_name in sample_names) {
    sample_rows <- final_data[final_data$Sample == sample_name, ]
    mat <- as.matrix(sample_rows[, c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")])
    mat <- log(mat)
    rownames(mat) <- paste0("event_", seq_len(nrow(mat)))
    ff <- flowCore::flowFrame(exprs = mat)
    fcs_path <- file.path(output_dir, paste0(sample_name, ".fcs"))
    flowCore::write.FCS(ff, filename = fcs_path)

    sample_meta_list[[length(sample_meta_list) + 1]] <- data.frame(
      Sample = sample_name,
      Treatment = sample_rows$Treatment[1],
      Batch = sample_rows$Batch[1],
      Setting = sample_rows$Setting[1],
      fcs_path = fcs_path,
      stringsAsFactors = FALSE
    )
  }

  sample_meta <- do.call(rbind, sample_meta_list)

  cell_meta <- final_data[, c("Sample", "Treatment", "Batch", "Setting", "CellType")]

  list(sample_meta = sample_meta, cell_meta = cell_meta)
}


#' Simulate differential expression (DE) flow cytometry data
#'
#' Generates simulated flow cytometry data with differential expression
#' of Marker2 in a target cell type across treatment groups, batch effects,
#' and multiple shift magnitudes. Writes per-sample FCS files and returns
#' sample- and cell-level metadata.
#'
#' @param groups Character vector of treatment group names.
#' @param batches Character vector of batch names.
#' @param sd_shifts Numeric vector of standard deviation shifts applied
#'   to Marker2 in the non-Baseline group for target cells.
#' @param samples_per_group Integer. Number of samples per group per shift.
#' @param mean_cells Integer. Mean number of cells per sample.
#' @param sd_cells Numeric. Standard deviation of cell count.
#' @param seed Integer. Random seed for reproducibility.
#' @param output_dir Character. Directory path where FCS files are written.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{sample_meta}{A \code{data.frame} with one row per sample
#'       and columns \code{Sample}, \code{Treatment}, \code{Batch},
#'       \code{SD_Shift}, \code{fcs_path}.}
#'     \item{cell_meta}{A \code{data.frame} with one row per cell
#'       and columns \code{Sample}, \code{Treatment}, \code{Batch},
#'       \code{SD_Shift}, \code{CellType}.}
#'   }
#'
#' @examples
#' \dontrun{
#' sim <- simulate_DE_data(
#'   samples_per_group = 2,
#'   mean_cells = 1000,
#'   sd_cells = 100,
#'   output_dir = file.path(tempdir(), "sim_DE_fcs")
#' )
#' head(sim$sample_meta)
#' head(sim$cell_meta)
#' }
#'
#' @export
simulate_DE_data <- function(groups = c("Baseline", "Activation"),
                             batches = c("Batch1", "Batch2"),
                             sd_shifts = c(0.5, 1, 2),
                             samples_per_group = 6,
                             mean_cells = 50000,
                             sd_cells = 500,
                             seed = 42,
                             output_dir = file.path(tempdir(), "sim_DE_fcs")) {

  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("Package 'flowCore' is required for simulate_DE_data(). ",
         "Install it with: BiocManager::install('flowCore')", call. = FALSE)
  }

  set.seed(seed)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize storage
  data_list <- list()

  # Simulate data
  for (sd_shift in sd_shifts) {
    for (group in groups) {
      for (sample_id in 1:samples_per_group) {
        batch <-
          if (sample_id %% 2 == 0) {
            "Batch2"
          } else {
            "Batch1"
          }
        total_cells <- round(stats::rnorm(1, mean = mean_cells, sd = sd_cells))
        total_cells <- max(total_cells, 1000)  # Ensure a minimum number of cells

        # Fixed proportion of cell type of interest
        proportion <- 0.05
        num_interest <- round(total_cells * proportion)
        num_other <- total_cells - num_interest
        cell_types <- c(rep("target", num_interest), rep("other", num_other))
        cell_types <- sample(cell_types)

        # Simulate expression data
        marker1 <- numeric(total_cells)
        marker2 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 1.5)
        marker3 <- stats::rlnorm(total_cells, meanlog = 0, sdlog = 2.5)
        marker4 <- numeric(total_cells)
        marker5 <- numeric(total_cells)

        # Assign Marker1, Marker4 Marker5 based on cell type
        marker1[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 2)
        marker1[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 5, sdlog = 2)  # Shift by 5 SD
        marker4[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.2)
        marker4[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 3, sdlog = 1.2)  # Shift by 3 SD
        marker5[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.8)
        marker5[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + 7, sdlog = 1.8)  # Shift by 7 SD

        # Marker2: differential expression between groups in activated only
        marker2[cell_types == "other"] <- stats::rlnorm(sum(cell_types == "other"), meanlog = 0, sdlog = 1.5)
        if (group == "Baseline") {
          marker2[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0, sdlog = 1.5)
        } else {
          marker2[cell_types == "target"] <- stats::rlnorm(sum(cell_types == "target"), meanlog = 0 + sd_shift, sdlog = 1.5)
        }

        # Add batch effect
        if (batch == "Batch2") {
          marker1 <- marker1 * stats::rlnorm(total_cells, meanlog = 0.1, sdlog = 0.5)
          marker2 <- marker2 * stats::rlnorm(total_cells, meanlog = 0.2, sdlog = 0.3)
          marker3 <- marker3 * stats::rlnorm(total_cells, meanlog = 0.3, sdlog = 0.4)
          marker4 <- marker4 * stats::rlnorm(total_cells, meanlog = 0.4, sdlog = 0.3)
          marker5 <- marker5 * stats::rlnorm(total_cells, meanlog = 0.5, sdlog = 0.5)
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
        data_list[[length(data_list) + 1]] <- df
      }
    }
  }

  # Combine
  final_data <- do.call(rbind, data_list) |>
    (\(x) {
      x$Treatment <- factor(x = x$Treatment,
                             levels = c("Baseline", "Activation"))
      x
    })()

  # Write FCS files and build sample_meta
  sample_names <- unique(final_data$Sample)
  sample_meta_list <- list()

  for (sample_name in sample_names) {
    sample_rows <- final_data[final_data$Sample == sample_name, ]
    mat <- as.matrix(sample_rows[, c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5")])
    mat <- log(mat)
    rownames(mat) <- paste0("event_", seq_len(nrow(mat)))
    ff <- flowCore::flowFrame(exprs = mat)
    fcs_path <- file.path(output_dir, paste0(sample_name, ".fcs"))
    flowCore::write.FCS(ff, filename = fcs_path)

    sample_meta_list[[length(sample_meta_list) + 1]] <- data.frame(
      Sample = sample_name,
      Treatment = as.character(sample_rows$Treatment[1]),
      Batch = sample_rows$Batch[1],
      SD_Shift = sample_rows$SD_Shift[1],
      fcs_path = fcs_path,
      stringsAsFactors = FALSE
    )
  }

  sample_meta <- do.call(rbind, sample_meta_list)

  cell_meta <- final_data[, c("Sample", "Treatment", "Batch", "SD_Shift", "CellType")]

  list(sample_meta = sample_meta, cell_meta = cell_meta)
}
