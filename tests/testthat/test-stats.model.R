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
  .tdr.obj <- list(
    cells = list(sample1 = "path1", sample2 = "path2"),
    map = list(fdens = matrix(runif(20), nrow = 10, ncol = 2))
  )
  .design <- matrix(1, nrow = 3, ncol = 1)  # Wrong number of rows
  
  expect_error(get.lm(.tdr.obj = .tdr.obj, .design = .design),
               "Number of rows in design matrix must be equal to the number of samples")
})

test_that("get.lm requires map", {
  .tdr.obj <- list(cells = list(sample1 = "path1"))
  .design <- matrix(1, nrow = 1, ncol = 1)
  
  expect_error(get.lm(.tdr.obj = .tdr.obj, .design = .design),
               "First run get.map")
})

# Test for get.pbDE

test_that("get.pbDE validates input correctly", {
  # Test error when design matrix doesn't match cells
  .tdr.obj <- list(
    cells = list(sample1 = "path1", sample2 = "path2"),
    assay.type = "RNA"
  )
  .design <- matrix(1, nrow = 3, ncol = 1)  # Wrong number of rows
  
  expect_error(get.pbDE(.tdr.obj = .tdr.obj, .design = .design),
               "Number of rows in design matrix must be equal to the number of samples")
})

test_that("get.pbDE validates geneset.ls for RNA", {
  .tdr.obj <- list(
    cells = list(sample1 = "path1"),
    assay.type = "cyto"
  )
  .design <- matrix(1, nrow = 1, ncol = 1)
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  expect_error(get.pbDE(.tdr.obj = .tdr.obj, .design = .design, .geneset.ls = .geneset.ls),
               ".geneset.ls is only supported for RNA assay type")
})

# Test for get.markerDE

test_that("get.markerDE validates geneset.ls correctly", {
  .tdr.obj <- list(
    assay.type = "cyto",
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .geneset.ls <- list(gene1 = c("A", "B"))
  
  expect_error(get.markerDE(.tdr.obj = .tdr.obj, .geneset.ls = .geneset.ls, .id1 = "A", .id2 = "B"),
               ".geneset.ls is only supported for RNA assay type")
})

test_that("get.markerDE requires valid id parameters", {
  .tdr.obj <- list(
    cells = list(sample1 = "path1"),
    assay.type = "RNA",
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  .design <- matrix(1, nrow = 1, ncol = 1)
  
  expect_error(get.markerDE(.tdr.obj = .tdr.obj),
               "Please provide .id1 or .id1.idx")
  expect_error(get.markerDE(.tdr.obj = .tdr.obj, .id1 = "C"),
               "not found in clustering")
})

test_that("get.marker is deprecated", {
  .tdr.obj <- list(
    assay.type = "cyto",
    graph = list(clustering = list(ids = factor(c("A", "B"))))
  )
  
  expect_warning(
    try(get.marker(.tdr.obj = .tdr.obj, .id1 = "A"), silent = TRUE),
    "deprecated"
  )
})
  
