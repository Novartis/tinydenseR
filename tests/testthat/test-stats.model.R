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
    map = list(norm = matrix(runif(20), nrow = 10, ncol = 2)),
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

# =============================================================================
# Mode detection tests
# =============================================================================

test_that("get.pbDE auto-detects design mode when .design is provided", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1", sample2 = "path2"),
    config = list(assay.type = "RNA")
  )
  .design <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2)  # rank-deficient
  # Should auto-detect design mode; error will be about design rank, not mode
  expect_error(get.pbDE(.tdr.obj, .design = .design),
               "not of full rank|not estimable|Number of rows")
})

test_that("get.pbDE auto-detects marker mode when .id is provided without .design", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  # Should auto-detect marker mode; will error later in the pipeline
  expect_error(
    get.pbDE(.tdr.obj, .id = "A"),
    "(?!Cannot determine analysis mode)",
    perl = TRUE
  )
})

test_that("get.pbDE errors when neither .design nor .id is provided", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA")
  )
  expect_error(get.pbDE(.tdr.obj),
               "Cannot determine analysis mode")
})

test_that("get.pbDE errors when .design AND .id2 are both provided", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1", sample2 = "path2"),
    config = list(assay.type = "RNA")
  )
  .design <- matrix(1, nrow = 2, ncol = 1)
  expect_error(
    get.pbDE(.tdr.obj, .design = .design, .id = "A", .id2 = "B"),
    "Conflicting arguments"
  )
})

test_that("get.pbDE errors in design mode without .design", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA")
  )
  expect_error(
    get.pbDE(.tdr.obj, .mode = "design"),
    "In design mode, .design is required"
  )
})

test_that("get.pbDE errors in marker mode with .design", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA")
  )
  .design <- matrix(1, nrow = 1, ncol = 1)
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker", .design = .design, .id = "A"),
    ".design is not used in marker mode"
  )
})

test_that("get.pbDE errors in marker mode without .id or .id.idx", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA")
  )
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker"),
    "In marker mode, .id or .id.idx"
  )
})

test_that("get.pbDE warns when design-mode-only args used in marker mode", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  expect_warning(
    try(get.pbDE(.tdr.obj, .mode = "marker", .id = "A",
                 .contrasts = matrix(1)), silent = TRUE),
    ".contrasts is not used in marker mode"
  )
  expect_warning(
    try(get.pbDE(.tdr.obj, .mode = "marker", .id = "A",
                 .block = "Donor"), silent = TRUE),
    ".block is not used in marker mode"
  )
})

# =============================================================================
# Marker mode validation (via get.pbDE)
# =============================================================================

test_that("get.pbDE marker mode validates geneset.ls correctly", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  expect_error(get.pbDE(.tdr.obj, .mode = "marker", .id = "A", .id2 = "B",
                        .geneset.ls = .geneset.ls),
               ".geneset.ls is only supported for RNA assay type")
})

test_that("get.pbDE marker mode validates .id labels", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_error(get.pbDE(.tdr.obj, .mode = "marker", .id = "C"),
               "not found in clustering")
})

# =============================================================================
# Deprecated alias tests (.population.name / .comparison.name)
# =============================================================================

test_that("get.pbDE warns when .population.name is used", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1", sample2 = "path2"),
    config = list(assay.type = "RNA")
  )
  .design <- matrix(1, nrow = 2, ncol = 1)
  expect_warning(
    try(get.pbDE(.tdr.obj, .design = .design,
                 .population.name = "test"), silent = TRUE),
    "`.population.name` is deprecated"
  )
})

test_that("get.pbDE warns when .comparison.name is used", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  expect_warning(
    try(get.pbDE(.tdr.obj, .mode = "marker", .id = "A",
                 .comparison.name = "test"), silent = TRUE),
    "`.comparison.name` is deprecated"
  )
})

# =============================================================================
# Backward-compatible get.markerDE tests (via deprecated wrapper)
# =============================================================================

test_that("get.markerDE is deprecated and dispatches to get.pbDE", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_warning(
    try(get.markerDE(.tdr.obj, .id1 = "A"), silent = TRUE),
    "soft-deprecated"
  )
})

test_that("get.markerDE validates geneset.ls correctly (backward compat)", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  # The deprecation warning and error should both fire
  expect_error(
    suppressWarnings(
      get.markerDE(.tdr.obj, .geneset.ls = .geneset.ls, .id1 = "A", .id2 = "B")
    ),
    ".geneset.ls is only supported for RNA assay type"
  )
})

test_that("get.markerDE requires valid id parameters (backward compat)", {
  .tdr.obj <- TDRObj(
    cells = list(sample1 = "path1"),
    config = list(assay.type = "RNA"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # .id1 = NULL should trigger "In marker mode, .id or .id.idx"
  expect_error(
    suppressWarnings(get.markerDE(.tdr.obj)),
    ".id or .id.idx"
  )
  # Invalid .id1 should trigger "not found in clustering"
  expect_error(
    suppressWarnings(get.markerDE(.tdr.obj, .id1 = "C")),
    "not found in clustering"
  )
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

# =============================================================================
# Regression: .id2.idx priority over .id2 parameter (via get.pbDE marker mode)
# =============================================================================

test_that("get.pbDE marker mode .id2.idx priority: invalid .id2 bypassed when .id2.idx provided", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Provide: valid .id, INVALID .id2 (not in clustering), valid .id2.idx
  # If priority is correct: should pass .id2 validation (uses .id2.idx instead)
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker",
             .id = "A", .id2 = "INVALID_CELLTYPE", .id2.idx = 1:5),
    regex = "(?!not found in clustering)",
    perl = TRUE
  )
})

test_that("get.pbDE marker mode .id2.idx consistency check", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Test 1: With .id2.idx, invalid .id2 should be bypassed
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker",
             .id = "A", .id2 = "INVALID", .id2.idx = 1:5),
    "(?!not found in clustering)",
    perl = TRUE
  )
  
  # Test 2: Without .id2.idx, invalid .id2 should be caught
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker",
             .id = "A", .id2 = "INVALID"),
    "INVALID.*not found in clustering"
  )
})

test_that("get.pbDE marker mode .id2.idx code path: priority at top level", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  # Provide both .id2 (special mode) and .id2.idx
  # .id2.idx should take priority, so "..all.other.landmarks.." message should not appear
  expect_error(
    get.pbDE(.tdr.obj, .mode = "marker",
             .id = "A", .id2 = "..all.other.landmarks..", .id2.idx = 1:5),
    "(?!using.*all\\.other\\.landmarks)",
    perl = TRUE
  )
})

# Also retain backward-compat versions of the same regression tests via get.markerDE
test_that("get.markerDE .id2.idx priority: invalid .id2 bypassed when .id2.idx provided", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_error(
    suppressWarnings(
      get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID_CELLTYPE", .id2.idx = 1:5)
    ),
    regex = "(?!not found in clustering)",
    perl = TRUE
  )
})

test_that("get.markerDE .id2.idx priority: API consistency check", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_error(
    suppressWarnings(
      get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID", .id2.idx = 1:5)
    ),
    "(?!not found in clustering)",
    perl = TRUE
  )
  
  expect_error(
    suppressWarnings(
      get.markerDE(.tdr.obj, .id1 = "A", .id2 = "INVALID")
    ),
    "INVALID.*not found in clustering"
  )
})

test_that("get.markerDE .id2.idx code path: priority at top level", {
  .tdr.obj <- TDRObj(
    config = list(assay.type = "cyto"),
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_error(
    suppressWarnings(
      get.markerDE(.tdr.obj, .id1 = "A", .id2 = "..all.other.landmarks..", .id2.idx = 1:5)
    ),
    "(?!using.*all\\.other\\.landmarks)",
    perl = TRUE
  )
})

