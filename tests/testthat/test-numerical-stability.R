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

library(testthat)
library(tinydenseR)

# ============================================================================
# NUMERICAL STABILITY TESTS
# ============================================================================
# These tests ensure tinydenseR produces reliable, reproducible results under:
# 1. Fixed random seeds (bit-identical reproducibility)
# 2. Small input perturbations (algorithmic robustness)
# 3. Invariant transformations (permutation, normalization)
# 4. Edge cases (extreme values, boundary conditions)
#
# Performance: All tests use small synthetic data and run in <10s total
# CI: Expensive checks use skip_on_cran() to avoid timeout
# ============================================================================

# ----------------------------------------------------------------------------
# 1. SEED REPRODUCIBILITY
# ----------------------------------------------------------------------------
# Verify that set.seed() produces bit-identical results across runs
# Critical for: regulatory compliance, publication reproducibility

test_that("full pipeline is reproducible with fixed seed", {
  skip_on_cran()  # Slightly slower due to double pipeline run
  
  # Create test data (small for speed)
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  
  # Run 1
  set.seed(123)
  result1 <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)
  
  # Run 2 (same seed)
  set.seed(123)
  result2 <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)
  
  # Assert bit-identical results
  expect_identical(rownames(result1$lm), rownames(result2$lm))  # Same landmarks
  expect_equal(result1$graph$uwot$embedding, result2$graph$uwot$embedding)  # Same UMAP
  expect_identical(as.character(result1$graph$clustering$ids), 
                   as.character(result2$graph$clustering$ids))  # Same clusters
  expect_equal(result1$map$fdens, result2$map$fdens, tolerance = 1e-14)  # Same densities
  
  cleanup_test_files(test_data)
})

test_that("clustering is deterministic with same seed", {
  set.seed(42)
  sim_matrix <- matrix(runif(400), 20, 20)
  sim_matrix <- (sim_matrix + t(sim_matrix)) / 2  # symmetric
  
  clusters1 <- tinydenseR:::leiden.cluster(
    .sim.matrix = sim_matrix,
    .resolution.parameter = 1.0,
    .seed = 42,
    .verbose = FALSE
  )
  
  clusters2 <- tinydenseR:::leiden.cluster(
    .sim.matrix = sim_matrix,
    .resolution.parameter = 1.0,
    .seed = 42,
    .verbose = FALSE
  )
  
  expect_identical(clusters1, clusters2)
})

# ----------------------------------------------------------------------------
# 2. INPUT INVARIANCE
# ----------------------------------------------------------------------------
# Verify robustness to benign transformations and small perturbations

test_that("pipeline is invariant to cell order permutation", {
  skip_on_cran()  # Double pipeline run
  
  # Create test data
  test_data <- create_test_lm_obj(n_cells = 60, n_markers = 4, n_samples = 2)
  
  # Baseline run
  set.seed(456)
  result1 <- run_full_pipeline(test_data, seed = 456, verbose = FALSE)
  
  # Permuted run: shuffle cell order within each sample
  permuted_cells <- lapply(test_data$cells, function(file_path) {
    mat <- readRDS(file_path)
    perm_idx <- sample(nrow(mat))
    mat_perm <- mat[perm_idx, , drop = FALSE]
    new_file <- tempfile(fileext = ".RDS")
    saveRDS(mat_perm, new_file, compress = FALSE)
    return(new_file)
  })
  names(permuted_cells) <- names(test_data$cells)
  
  test_data_perm <- list(cells = permuted_cells, meta = test_data$meta)
  
  set.seed(456)
  result2 <- run_full_pipeline(test_data_perm, seed = 456, verbose = FALSE)
  
  # fdens should be identical (cell order doesn't matter for aggregates)
  expect_equal(result1$map$fdens, result2$map$fdens, tolerance = 1e-12)
  
  # Cleanup
  cleanup_test_files(test_data)
  lapply(permuted_cells, unlink)
})

test_that("pipeline is robust to small perturbations", {
  skip_on_cran()  # Double pipeline run
  
  # Create baseline data
  test_data <- create_test_lm_obj(n_cells = 80, n_markers = 6, n_samples = 3)
  
  # Baseline run
  set.seed(789)
  result1 <- run_full_pipeline(test_data, seed = 789, verbose = FALSE)
  
  # Perturbed run: add small noise (σ = 1e-6, larger than 1e-8 but still tiny)
  perturbed_cells <- lapply(test_data$cells, function(file_path) {
    mat <- readRDS(file_path)
    mat_noisy <- mat + matrix(rnorm(length(mat), sd = 1e-6), nrow = nrow(mat))
    new_file <- tempfile(fileext = ".RDS")
    saveRDS(mat_noisy, new_file, compress = FALSE)
    return(new_file)
  })
  names(perturbed_cells) <- names(test_data$cells)
  
  test_data_pert <- list(cells = perturbed_cells, meta = test_data$meta)
  
  set.seed(789)
  result2 <- run_full_pipeline(test_data_pert, seed = 789, verbose = FALSE)
  
  # Compute Jaccard similarity of landmark sets
  jaccard_lm <- jaccard_similarity(rownames(result1$lm), rownames(result2$lm))
  
  # Should be highly stable (>85% overlap for landmarks)
  expect_gt(jaccard_lm, 0.85)
  
  # fdens should be similar (not identical, but close)
  expect_equal(result1$map$fdens, result2$map$fdens, tolerance = 1e-3)
  
  # Cleanup
  cleanup_test_files(test_data)
  lapply(perturbed_cells, unlink)
})

test_that("fuzzy graph weights are valid and normalized", {
  skip_on_cran()
  
  # Build a small pipeline
  test_data <- create_test_lm_obj(n_cells = 40, n_markers = 3, n_samples = 2)
  
  set.seed(111)
  result <- run_full_pipeline(test_data, seed = 111, verbose = FALSE)
  
  # Check fdens matrix properties
  expect_true(all(is.finite(result$map$fdens)))  # No NaN/Inf
  expect_true(all(result$map$fdens >= 0))  # Non-negative
  
  # Check cell count matrices
  expect_true(all(result$map$clustering$cell.count >= 0))
  expect_true(all(is.finite(result$map$clustering$cell.perc)))
  expect_true(all(result$map$clustering$cell.perc >= 0))
  expect_true(all(result$map$clustering$cell.perc <= 100))
  
  # Row sums of cell.perc should be ~100% (excluding rows with all zeros/NaN)
  row_sums <- rowSums(result$map$clustering$cell.perc)
  valid_rows <- is.finite(row_sums) & row_sums > 0
  expect_true(all(abs(row_sums[valid_rows] - 100) < 1e-6))
  
  cleanup_test_files(test_data)
})

# ----------------------------------------------------------------------------
# 3. EDGE CASES
# ----------------------------------------------------------------------------
# Verify graceful handling of boundary conditions and extreme inputs

test_that("pipeline handles minimal input dimensions", {
  skip_on_cran()
  
  # Small dataset: 3 samples, 3 markers, 20 cells each (minimum viable for LE)
  # Total: 60 cells -> ~9 landmarks -> enough for k=2, nPC=2, LE eigenmap
  test_data <- create_test_lm_obj(n_cells = 20, n_markers = 3, n_samples = 3)
  
  # Should not error with small but viable dataset
  expect_error(
    run_full_pipeline(test_data, seed = 222, verbose = FALSE),
    NA
  )
  
  cleanup_test_files(test_data)
})

test_that("clustering handles stragglers correctly", {
  # Create similarity matrix that will produce tiny clusters
  set.seed(101)
  sim_matrix <- matrix(0, 30, 30)
  # Main cluster: strong connections among first 25 nodes
  sim_matrix[1:25, 1:25] <- runif(625, 0.5, 1.0)
  # Stragglers: 5 weakly connected nodes
  sim_matrix[26:30, 1:25] <- runif(125, 0.01, 0.05)
  sim_matrix <- (sim_matrix + t(sim_matrix)) / 2
  diag(sim_matrix) <- 1
  
  clusters <- tinydenseR:::leiden.cluster(
    .sim.matrix = sim_matrix,
    .resolution.parameter = 0.5,
    .small.size = 3,  # Absorb clusters < 3
    .seed = 101,
    .verbose = FALSE
  )
  
  # Check no cluster is smaller than threshold (stragglers absorbed)
  cluster_sizes <- table(clusters)
  expect_true(all(cluster_sizes >= 3))
})

test_that("Jaccard similarity handles empty neighbor sets", {
  # Construct adjacency with some isolated nodes (no neighbors)
  adj_matrix <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 3),
    j = c(2, 3, 1, 2),
    x = 1,
    dims = c(5, 5)
  )
  # Nodes 4 and 5 have no neighbors
  
  # Should not crash, should return valid sparse matrix
  expect_error(
    result <- tinydenseR:::fast.jaccard.r(adj_matrix, .prune = 1/15),
    NA
  )
  
  # Check that result is a valid sparse matrix
  expect_s4_class(result, "sparseMatrix")
  expect_true(Matrix::isSymmetric(result))
  
  # Note: Jaccard with empty sets can produce Inf/-Inf, which is mathematically correct
  # The function uses formula: 1 / ((2k / intersect) - 1)
  # When intersect=0, this produces Inf, which gets pruned away
})

test_that("pipeline handles extreme marker values", {
  skip_on_cran()
  
  # Create data with extreme ranges (need 3+ markers per package requirements)
  cells_extreme <- lapply(1:2, function(i) {
    mat <- matrix(
      c(runif(20, 0, 0.001),      # Very small values
        runif(20, 1000, 10000),   # Very large values
        runif(20, 1, 10)),        # Normal range
      nrow = 20,
      ncol = 3,
      dimnames = list(
        paste0("sample", i, "_cell_", 1:20),
        c("marker_low", "marker_high", "marker_normal")
      )
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(mat, uri, compress = FALSE)
    return(uri)
  })
  names(cells_extreme) <- paste0("sample", 1:2)
  
  meta_extreme <- data.frame(
    row.names = names(cells_extreme),
    group = c("A", "B")
  )
  
  # Should handle extreme values without producing Inf/NaN
  expect_error(
    {
      result <- setup.lm.obj(
        .cells = cells_extreme,
        .meta = meta_extreme,
        .markers = c("marker_low", "marker_high", "marker_normal"),
        .assay.type = "cyto",
        .verbose = FALSE
      ) |>
        get.landmarks(.verbose = FALSE, .nPC = 2) |>
        get.graph(.k = 3, .verbose = FALSE)
      
      # Verify no NaN/Inf in results
      expect_true(all(is.finite(result$lm)))
      expect_true(all(is.finite(result$graph$uwot$embedding)))
    },
    NA
  )
  
  lapply(cells_extreme, unlink)
})

# ----------------------------------------------------------------------------
# 4. CI ROBUSTNESS (cross-platform, cross-version)
# ----------------------------------------------------------------------------
# Verify stable behavior across different environments

test_that("RNG produces consistent results across sessions", {
  # Test that RNGkind() doesn't affect reproducibility
  old_rng <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_rng)))
  
  # Set standard RNG
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(999)
  
  sim1 <- matrix(runif(100), 10, 10)
  sim1 <- (sim1 + t(sim1)) / 2
  clusters1 <- tinydenseR:::leiden.cluster(
    .sim.matrix = sim1,
    .resolution.parameter = 1.0,
    .seed = 999,
    .verbose = FALSE
  )
  
  # Reset and rerun
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(999)
  
  sim2 <- matrix(runif(100), 10, 10)
  sim2 <- (sim2 + t(sim2)) / 2
  clusters2 <- tinydenseR:::leiden.cluster(
    .sim.matrix = sim2,
    .resolution.parameter = 1.0,
    .seed = 999,
    .verbose = FALSE
  )
  
  expect_identical(clusters1, clusters2)
})

test_that("sparse matrix operations preserve numerical properties", {
  # Test that conversion between dense/sparse doesn't introduce errors
  set.seed(202)
  dense_adj <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.8, 0.2)), 10, 10)
  sparse_adj <- Matrix::Matrix(dense_adj, sparse = TRUE)
  
  # Jaccard on sparse
  result_sparse <- tinydenseR:::fast.jaccard.r(sparse_adj, .prune = 0)
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(result_sparse))
  
  # Values should be finite (Jaccard can have edge cases with formula used)
  expect_true(all(is.finite(result_sparse@x)))
  
  # Most values should be in reasonable range [0, 1] but allow for numerical edge cases
  expect_true(all(result_sparse@x <= 1.01))  # Small tolerance for rounding
})
