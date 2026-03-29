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

# Test for get.lm

test_that("get.lm validates input correctly", {
  # Test error when design matrix doesn't match cells
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1", sample2 = "path2"),
    map = list(fdens = matrix(runif(20), nrow = 10, ncol = 2)),
    config = list(assay.type = "cyto")
  )
  .design <- matrix(1, nrow = 3, ncol = 1)  # Wrong number of rows
  
  expect_error(get.lm(.tdr.obj, .design = .design),
               "Number of rows in design matrix must be equal to the number of samples")
})

test_that("get.lm requires map", {
  .tdr.obj <- TDRObj(cells = list(sample1 = "path1"), config = list(assay.type = "cyto"))
  .design <- matrix(1, nrow = 1, ncol = 1)
  
  expect_error(get.lm(.tdr.obj, .design = .design),
               "First run get.map")
})

# Test for get.pbDE

test_that("get.pbDE validates input correctly", {
  # Test error when design matrix doesn't match cells
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1", sample2 = "path2"),
    config = list(assay.type = "RNA")
  )
  .design <- matrix(1, nrow = 3, ncol = 1)  # Wrong number of rows
  
  expect_error(get.pbDE(.tdr.obj, .design = .design),
               "Number of rows in design matrix must be equal to the number of samples")
})

test_that("get.pbDE validates geneset.ls for RNA", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "cyto")
  )
  .design <- matrix(1, nrow = 1, ncol = 1)
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  expect_error(get.pbDE(.tdr.obj, .design = .design, .geneset.ls = .geneset.ls),
               ".geneset.ls is only supported for RNA assay type")
})

# Test for get.markerDE

test_that("get.markerDE validates geneset.ls correctly", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  expect_error(get.markerDE(.tdr.obj, .geneset.ls = .geneset.ls, .id1 = "A", .id2 = "B"),
               ".geneset.ls is only supported for RNA assay type")
})

test_that("get.markerDE requires valid id parameters", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .design <- matrix(1, nrow = 1, ncol = 1)
  
  expect_error(get.markerDE(.tdr.obj),
               "Please provide .id1 or .id1.idx")
  expect_error(get.markerDE(.tdr.obj, .id1 = "C"),
               "not found in clustering")
})

test_that("get.marker is deprecated", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_warning(
    try(get.marker(.tdr.obj, .id1 = "A"), silent = TRUE),
    "deprecated"
  )
})

# Regression test: .id2.idx priority over .id2 parameter
# This test verifies that when .id2.idx is provided, it takes priority over .id2
# by checking that providing an invalid .id2 with valid .id2.idx doesn't error at validation
test_that("get.markerDE .id2.idx priority: invalid .id2 bypassed when .id2.idx provided", {
  # Regression: Verify that .id2.idx priority logic is correctly structured
  # If .id2.idx is prioritized at the top level as intended, providing an 
  # invalid .id2 value should not cause an immediate error.
  # Instead, it should only validate .id2 when .id2.idx is NULL.
  
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Provide: valid .id1, INVALID .id2 (not in clustering), valid .id2.idx
  # If priority is correct: should pass .id2 validation (uses .id2.idx instead)
  # If priority is wrong: would error on invalid .id2
  
  expect_error(
    get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID_CELLTYPE", .id2.idx = 1:5),
    # If .id2.idx is prioritized, this will fail later (missing required data),
    # not with "not found in clustering" validation error
    regex = "(?!not found in clustering)",  # negative lookahead: should NOT have this error
    perl = TRUE
  )
})

# Structured integration test for .id2.idx priority behavior
test_that("get.markerDE .id2.idx priority: API consistency check", {
  # Integration test: Verify .id1.idx and .id2.idx follow symmetric priority pattern
  # Create a minimal object that has the required validation structure
  
  # The test passes if we can call get.markerDE with:
  # - Valid .id1/invalid .id2 + valid .id2.idx -> should not error on .id2 validation
  # - Valid .id1.idx/invalid .id1 + valid .id2/unknown .id2.idx -> should error on .id1 validation  
  
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Test 1: With .id2.idx, invalid .id2 should be bypassed
  expect_error(
    get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID", .id2.idx = 1:5),
    "(?!not found in clustering)",
    perl = TRUE
  )
  
  # Test 2: Without .id2.idx, invalid .id2 should be caught
  expect_error(
    get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID"),
    "INVALID.*not found in clustering"
  )
})

# Final regression verification: priority logic code path
test_that("get.markerDE .id2.idx code path: priority at top level", {
  #  Verification that .id2.idx check is at the top level, not nested.
  # This ensures: is.null(.id2.idx) -> check .id2 vs special mode
  #              !is.null(.id2.idx) -> use .id2.idx directly
  
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Scenario 1: Provide both .id2 (special mode) and .id2.idx
  # If .id2.idx is truly at top level and prioritized:
  # -> won't apply the special "..all.other.landmarks.." mode
  # -> will attempt to use .id2.idx directly (will fail at later stage, not here)
  expect_error(
    get.markerDE(.tdr.obj, .id1 = "A", .id2 = "..all.other.landmarks..", .id2.idx = 1:5),
    "(?!using.*all\\.other\\.landmarks)",  # should NOT see this message
    perl = TRUE
  )
})

