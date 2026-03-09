#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ======================================================================
# simulate_DA_data
# ======================================================================

test_that("simulate_DA_data returns correct structure", {
  skip_if_not_installed("flowCore")

  res <- simulate_DA_data(samples_per_group = 2,
                          mean_cells = 500,
                          sd_cells = 50,
                          seed = 42)
  withr::defer(unlink(dirname(res$sample_meta$fcs_path[1]),
                      recursive = TRUE))

  expect_type(res, "list")
  expect_named(res, c("sample_meta", "cell_meta"))

 expect_true(all(c("Sample", "Treatment", "Batch", "Setting", "fcs_path")
                  %in% colnames(res$sample_meta)))
  expect_true(all(c("Sample", "Treatment", "Batch", "Setting", "CellType")
                  %in% colnames(res$cell_meta)))

  # Default: 2 groups * 2 samples_per_group * 3 settings = 12
  expect_equal(nrow(res$sample_meta), 2L * 2L * 3L)

  expect_true(all(file.exists(res$sample_meta$fcs_path)))
})

test_that("simulate_DA_data FCS files are valid", {
  skip_if_not_installed("flowCore")

  res <- simulate_DA_data(samples_per_group = 2,
                          mean_cells = 500,
                          sd_cells = 50,
                          seed = 42)
  withr::defer(unlink(dirname(res$sample_meta$fcs_path[1]),
                      recursive = TRUE))

  ff <- flowCore::read.FCS(res$sample_meta$fcs_path[1])
  expr_mat <- flowCore::exprs(ff)

  expect_true(all(paste0("Marker", 1:5) %in% colnames(expr_mat)))
  expect_true(all(is.finite(expr_mat)))
})

test_that("simulate_DA_data is deterministic", {
  skip_if_not_installed("flowCore")

  res1 <- simulate_DA_data(samples_per_group = 2,
                           mean_cells = 500,
                           sd_cells = 50,
                           seed = 42)
  withr::defer(unlink(dirname(res1$sample_meta$fcs_path[1]),
                      recursive = TRUE))

  res2 <- simulate_DA_data(samples_per_group = 2,
                           mean_cells = 500,
                           sd_cells = 50,
                           seed = 42)
  withr::defer(unlink(dirname(res2$sample_meta$fcs_path[1]),
                      recursive = TRUE))

  expect_identical(res1$cell_meta, res2$cell_meta)

  ff1 <- flowCore::read.FCS(res1$sample_meta$fcs_path[1])
  ff2 <- flowCore::read.FCS(res2$sample_meta$fcs_path[1])
  expect_equal(flowCore::exprs(ff1), flowCore::exprs(ff2))
})

# ======================================================================
# simulate_DE_data
# ======================================================================

test_that("simulate_DE_data returns correct structure", {
  skip_if_not_installed("flowCore")

  res <- simulate_DE_data(samples_per_group = 2,
                          mean_cells = 500,
                          sd_cells = 50,
                          seed = 42)
  withr::defer(unlink(dirname(res$sample_meta$fcs_path[1]),
                      recursive = TRUE))

  expect_type(res, "list")
  expect_named(res, c("sample_meta", "cell_meta"))

  expect_true(all(c("Sample", "Treatment", "Batch", "SD_Shift", "fcs_path")
                  %in% colnames(res$sample_meta)))
  expect_true(all(c("Sample", "Treatment", "Batch", "SD_Shift", "CellType")
                  %in% colnames(res$cell_meta)))

  # Default: 2 groups * 2 samples_per_group * 3 sd_shifts = 12
  expect_equal(nrow(res$sample_meta), 2L * 2L * 3L)

  expect_true(all(file.exists(res$sample_meta$fcs_path)))
})
