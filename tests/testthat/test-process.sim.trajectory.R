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

test_that("fetch_trajectory_data works and returns expected structure", {
  skip_if_not(curl::has_internet(), "No internet connection")
  
  # Test the function
  trajectory_data <- fetch_trajectory_data()
  
  # Check structure
  expect_type(trajectory_data, "list")
  expect_named(trajectory_data, c("meta", "SCE"))
  
  # Check metadata
  expect_true(is.data.frame(trajectory_data$meta))
  expect_true(all(c("Condition", "Replicate", "Sample") %in% colnames(trajectory_data$meta)))
  
  # Check SCE object
  expect_s4_class(trajectory_data$SCE, "SingleCellExperiment")
  
  # Check dimensions make sense
  expect_gt(nrow(trajectory_data$SCE), 0)
  expect_gt(ncol(trajectory_data$SCE), 0)
  
  # Test that we can use it with our functions
  # Create .meta object containing sample-level data using get.meta
  .meta <- get.meta(.obj = trajectory_data$SCE,
                    .sample.var = "Sample",
                    .verbose = FALSE)
  
  cells <- get.cells.SCE(.sce.obj = trajectory_data$SCE, .meta = .meta, .sample.var = "Sample")
  expect_type(cells, "list")
  expect_length(cells, 6)  # Should have 6 samples (A_R1, A_R2, A_R3, B_R1, B_R2, B_R3)
})

test_that("fetch_trajectory_data fails gracefully without internet", {
  skip_if(curl::has_internet(), "Internet connection available")
  
  expect_error(fetch_trajectory_data(), "Internet connection required")
})
