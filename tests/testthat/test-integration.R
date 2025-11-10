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

# Integration tests for typical workflow

test_that("setup.lm.obj integration works", {
  # Test just the setup step which is more reliable
  test_data <- create_test_lm_obj(n_cells = 10, n_markers = 3, n_samples = 2)
  
  # Step 1: Setup
  .lm.obj <- setup.lm.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = c("marker_1", "marker_2", "marker_3"),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  
  expect_type(object = .lm.obj, type = "list")
  expect_true("cells" %in% names(x = .lm.obj))
  expect_true("metadata" %in% names(x = .lm.obj))
  expect_true("assay.type" %in% names(x = .lm.obj))
  expect_equal(.lm.obj$assay.type, "cyto")
  
  # Clean up temp files
  cleanup_test_files(test_data)
})

test_that("workflow functions validate inputs properly", {
  # Test that functions properly validate their inputs
  
  # Test setup.lm.obj validation with mismatched names
  cells_list <- list("sample1" = matrix(1:10, ncol=2))
  meta_df <- data.frame(condition = "A", row.names = "sample2")
  expect_error(setup.lm.obj(.cells = cells_list, .meta = meta_df),
               "Sample names mismatch between .cells and .meta")
  
  # Test get.landmarks validation with incomplete object
  .lm.obj <- list(assay.type = "RNA")
  expect_error(get.landmarks(.lm.obj = .lm.obj),
               "Invalid .lm.obj")  # Should error due to missing required fields
})
