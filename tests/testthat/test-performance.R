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

# Performance and boundary tests

test_that("functions handle empty inputs gracefully", {
  # Test with mismatched names (should trigger early validation)
  empty_cells <- list("sample1" = tempfile())
  empty_meta <- data.frame(condition = "A", row.names = "sample2")
  
  expect_error(setup.lm.obj(.cells = empty_cells, .meta = empty_meta),
               "names of .cells must be the same as rownames of .meta")
  
  # Test get.stats with invalid inputs
  expect_error(get.stats(.lm.obj = list()))
})

test_that("functions handle large datasets", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  
  # Test setup.lm.obj with larger mock data
  large_cells <- list(
    sample1 = tempfile(),
    sample2 = tempfile(),
    sample3 = tempfile()
  )
  
  # Create temporary files with larger matrices that have proper rownames
  for(i in 1:3) {
    mat <- matrix(runif(300), ncol=3)
    rownames(mat) <- paste0("cell_", 1:100)
    colnames(mat) <- paste0("marker_", 1:3)
    saveRDS(mat, file = large_cells[[i]])
  }
  
  large_meta <- data.frame(
    condition = c("A", "B", "C"),
    row.names = paste0("sample", 1:3)
  )
  
  result <- setup.lm.obj(
    .cells = large_cells, 
    .meta = large_meta, 
    .verbose = FALSE
  )
  
  expect_type(object = result, type = "list")
  expect_true(length(result$cells) == 3)
  
  # Clean up
  lapply(large_cells, unlink)
})

test_that("plotting functions handle different data sizes", {
  # Test with minimal data
  .lm.obj_small <- list(
    graph = list(
      uwot = list(embedding = matrix(runif(4), ncol=2)),
      clustering = list(ids = factor(c("A", "B")))
    )
  )
  
  result_small <- plotUMAP(.lm.obj = .lm.obj_small)
  expect_true("ggplot" %in% class(result_small))
  
  # Test with larger data
  .lm.obj_large <- list(
    graph = list(
      uwot = list(embedding = matrix(runif(2000), ncol=2)),
      clustering = list(ids = factor(sample(letters[1:5], 1000, replace = TRUE)))
    )
  )
  
  result_large <- plotUMAP(.lm.obj = .lm.obj_large)
  expect_true("ggplot" %in% class(result_large))
})

test_that("functions handle edge cases", {
  # Test with invalid data types
  expect_error(get.stats(.lm.obj = "not_a_list"))
  expect_error(get.dea(.lm.obj = NULL))
  
  # Test plotting functions with minimal valid data
  minimal_obj <- list(
    graph = list(
      uwot = list(embedding = matrix(c(1, 2, 3, 4), ncol=2)),
      clustering = list(ids = factor(c("A", "B")))
    )
  )
  
  result <- plotUMAP(.lm.obj = minimal_obj)
  expect_true("ggplot" %in% class(result))
})

test_that("functions validate input types properly", {
  # Test that functions check for proper input types
  expect_error(setup.lm.obj(.cells = "not_a_list", .meta = data.frame()))
  
  # Test with valid cells but wrong meta type - create a temporary RDS file with minimal data
  temp_file <- tempfile(fileext = ".RDS")
  minimal_matrix <- matrix(1:10, nrow = 2, ncol = 5, 
                           dimnames = list(c("cell1", "cell2"), c("gene1", "gene2", "gene3", "gene4", "gene5")))
  saveRDS(minimal_matrix, file = temp_file)
  empty_cells <- list("sample1" = temp_file)
  expect_error(setup.lm.obj(.cells = empty_cells, .meta = "not_a_dataframe"), 
               ".meta must be a data.frame object")
  unlink(temp_file)  # Clean up
  
  # Test with proper types but wrong structure
  wrong_cells <- list(1, 2, 3)  # numeric instead of file paths
  proper_meta <- data.frame(condition = "A", row.names = "sample1")
  expect_error(setup.lm.obj(.cells = wrong_cells, .meta = proper_meta))
})
