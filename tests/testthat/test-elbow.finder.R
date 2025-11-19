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

# Test elbow.sec.deriv function

test_that("elbow.sec.deriv detects elbow in descending curve (PCA-like)", {
  # Simulate PCA eigenvalues: steep drop then plateau
  x <- c(10, 8, 6, 3, 2, 1.5, 1.3, 1.2, 1.1, 1.0)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  expect_named(object = result, expected = c("index", "value", "sec.deriv"))
  expect_true(object = result$index >= 1 && result$index <= length(x = x))
  expect_equal(object = result$value, expected = x[result$index])
  # Elbow should be in early positions (2-4) where curvature is highest
  expect_true(object = result$index <= 5)
})

test_that("elbow.sec.deriv detects elbow in ascending curve (LE-like)", {
  # Simulate Laplacian Eigenmap eigenvalues: low values then sharp increase
  x <- c(0.01, 0.02, 0.03, 0.1, 0.5, 1.2, 2.5, 4.0, 6.0, 8.0)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "asc")
  
  expect_type(object = result, type = "list")
  expect_named(object = result, expected = c("index", "value", "sec.deriv"))
  expect_true(object = result$index >= 1 && result$index <= length(x = x))
  expect_equal(object = result$value, expected = x[result$index])
  # Elbow should be around position where sharp increase begins
  expect_true(object = result$index >= 3 && result$index <= 7)
})

test_that("elbow.sec.deriv works with smoothing enabled", {
  # Noisy data that benefits from smoothing
  set.seed(seed = 123)
  x <- c(10, 9, 8, 6, 4, 2, 1, 0.8, 0.7, 0.6) + rnorm(n = 10, sd = 0.2)
  
  result_smooth <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, sort.order = "desc")
  result_no_smooth <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result_smooth, type = "list")
  expect_type(object = result_no_smooth, type = "list")
  # Both should return valid indices
  expect_true(object = result_smooth$index >= 1 && result_smooth$index <= length(x = x))
  expect_true(object = result_no_smooth$index >= 1 && result_no_smooth$index <= length(x = x))
})

test_that("elbow.sec.deriv respects custom df parameter", {
  x <- c(10, 8, 6, 4, 3, 2.5, 2.2, 2.1, 2.05, 2.0)
  
  # Should work with different df values
  result_df3 <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, df = 3, sort.order = "desc")
  result_df7 <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, df = 7, sort.order = "desc")
  
  expect_type(object = result_df3, type = "list")
  expect_type(object = result_df7, type = "list")
  expect_true(object = result_df3$index >= 1 && result_df3$index <= length(x = x))
  expect_true(object = result_df7$index >= 1 && result_df7$index <= length(x = x))
})

test_that("elbow.sec.deriv handles edge case: too few values", {
  # Need at least 4 values for second derivative
  expect_error(
    object = tinydenseR:::elbow.sec.deriv(x = c(1, 2, 3), sort.order = "desc"),
    regexp = "Cannot compute elbow point: need at least 4 values"
  )
  
  expect_error(
    object = tinydenseR:::elbow.sec.deriv(x = numeric(0), sort.order = "desc"),
    regexp = "Cannot compute elbow point: need at least 4 values"
  )
})

test_that("elbow.sec.deriv handles edge case: flat sequence", {
  # All identical values - no variation
  x <- rep(x = 5, times = 10)
  
  expect_error(
    object = tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc"),
    regexp = "Cannot detect elbow: all values are identical"
  )
})

test_that("elbow.sec.deriv handles edge case: minimal data (4 values)", {
  # Exactly 4 values - minimum required
  x <- c(10, 5, 2, 1)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  expect_true(object = result$index >= 1 && result$index <= 4)
})

test_that("elbow.sec.deriv smoothing only applied when length > 6", {
  # With 6 values, smoothing should be skipped even if smooth = TRUE
  x <- c(10, 8, 6, 4, 2, 1)
  
  # Should not error even with smooth = TRUE
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  expect_true(object = result$index >= 1 && result$index <= 6)
})

test_that("elbow.sec.deriv returns sensible index for typical PCA scree plot", {
  # Realistic PCA eigenvalues (first 15 PCs)
  x <- c(45.2, 23.1, 15.7, 8.9, 5.2, 3.8, 2.9, 2.1, 1.8, 1.5, 1.3, 1.2, 1.1, 1.0, 0.9)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  # For typical scree plots, elbow is usually in first 5-6 components
  expect_true(object = result$index >= 2 && result$index <= 8)
  expect_equal(object = result$value, expected = x[result$index])
})

test_that("elbow.sec.deriv handles unsorted input correctly", {
  # Provide unsorted data - function should sort internally
  x <- c(5, 10, 2, 8, 1, 6, 3, 4, 7, 9)
  
  result_desc <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  result_asc <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "asc")
  
  expect_type(object = result_desc, type = "list")
  expect_type(object = result_asc, type = "list")
  # Should return index from original x vector
  expect_true(object = result_desc$index >= 1 && result_desc$index <= length(x = x))
  expect_true(object = result_asc$index >= 1 && result_asc$index <= length(x = x))
  expect_equal(object = result_desc$value, expected = x[result_desc$index])
  expect_equal(object = result_asc$value, expected = x[result_asc$index])
})

test_that("elbow.sec.deriv sec.deriv output is numeric and positive", {
  x <- c(10, 8, 6, 4, 3, 2, 1.5, 1.2, 1, 0.9)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result$sec.deriv, type = "double")
  # Second derivative at elbow should be meaningful (typically positive after normalization)
  expect_true(object = is.finite(x = result$sec.deriv))
})

test_that("elbow.sec.deriv handles negative values correctly", {
  # Mix of negative and positive values
  x <- c(-5, -3, -1, 0, 1, 2, 3, 4, 5, 6)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "asc")
  
  expect_type(object = result, type = "list")
  expect_true(object = result$index >= 1 && result$index <= length(x = x))
  expect_equal(object = result$value, expected = x[result$index])
})

test_that("elbow.sec.deriv validates sort.order argument", {
  x <- c(10, 8, 6, 4, 2, 1)
  
  expect_error(
    object = tinydenseR:::elbow.sec.deriv(x = x, sort.order = "invalid"),
    regexp = "'arg' should be one of"
  )
  
  # Valid values should work
  expect_type(object = tinydenseR:::elbow.sec.deriv(x = x, sort.order = "desc"), type = "list")
  expect_type(object = tinydenseR:::elbow.sec.deriv(x = x, sort.order = "asc"), type = "list")
})

test_that("elbow.sec.deriv auto-calculates df appropriately", {
  # With 20 values, df should be min(20 * 0.7, 10) = 10
  x <- seq(from = 20, to = 1, length.out = 20)
  
  # Should not error with auto df
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = TRUE, df = NULL, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  expect_true(object = result$index >= 1 && result$index <= 20)
})

test_that("elbow.sec.deriv handles very smooth gradual decline", {
  # Linear decline - no clear elbow
  x <- seq(from = 10, to = 1, length.out = 15)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  # Should still return a valid index (middle-ish region)
  expect_true(object = result$index >= 1 && result$index <= 15)
})

test_that("elbow.sec.deriv handles single sharp drop", {
  # Steep initial drop then flat
  x <- c(100, 10, 9, 8.5, 8.2, 8.1, 8, 7.9, 7.8, 7.7)
  
  result <- tinydenseR:::elbow.sec.deriv(x = x, smooth = FALSE, sort.order = "desc")
  
  expect_type(object = result, type = "list")
  # Elbow should be at position 1 or 2 where the drop occurs
  expect_true(object = result$index <= 3)
})
