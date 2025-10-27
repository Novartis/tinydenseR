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

# Simplified plotting tests focusing on input validation only

test_that("plotting functions validate input types", {
  # Test that plotting functions require proper input types
  expect_error(plotUMAP(.lm.obj = NULL))
  expect_error(plotUMAP(.lm.obj = "not_a_list"))
  expect_error(plotUMAP(.lm.obj = list()))
  
  expect_error(plotPCA(.lm.obj = NULL))
  expect_error(plotPCA(.lm.obj = "not_a_list"))
  
  expect_error(scatterPlot(.x.feature = NULL, .y.feature = NULL))
  # Test scatterPlot with mismatched vector lengths
  expect_error(scatterPlot(.x.feature = 1:5, .y.feature = 1:3))
})

test_that("plotUMAP works with minimal valid structure", {
  # Test plotUMAP with minimal required structure
  .lm.obj <- list(
    graph = list(
      uwot = list(embedding = matrix(c(1, 2, 3, 4), ncol=2)),
      clustering = list(ids = factor(c("A", "B")))
    )
  )
  
  result <- plotUMAP(.lm.obj = .lm.obj)
  expect_true("ggplot" %in% class(result))
})

test_that("scatterPlot works with numeric vectors", {
  # Test scatterPlot with simple numeric data
  .x.feature <- rnorm(10)
  .y.feature <- rnorm(10)
  
  result <- scatterPlot(.x.feature = .x.feature, .y.feature = .y.feature)
  expect_true("ggplot" %in% class(result))
})

test_that("plotting functions handle mismatched dimensions", {
  # Test scatterPlot with mismatched vector lengths
  expect_error(scatterPlot(.x.feature = 1:5, .y.feature = 1:3))
  
  # Test with zero-length vectors
  expect_error(scatterPlot(.x.feature = numeric(0), .y.feature = numeric(0)))
})

test_that("plotting functions validate required fields", {
  # Test plotUMAP missing required fields
  incomplete_obj1 <- list(graph = list())
  expect_error(plotUMAP(.lm.obj = incomplete_obj1))
  
  incomplete_obj2 <- list(graph = list(uwot = list()))
  expect_error(plotUMAP(.lm.obj = incomplete_obj2))
  
  incomplete_obj3 <- list(graph = list(uwot = list(embedding = matrix(c(1,2), ncol=2))))
  expect_error(plotUMAP(.lm.obj = incomplete_obj3))  # Missing clustering
})

test_that("input validation for other plotting functions", {
  # Test basic input validation for functions that require complex structures
  # Rather than mock the complex structures, just test they expect the right input types
  
  expect_error(plotPCA(.lm.obj = list()))
  expect_error(plotBeeswarm(.lm.obj = list()))
  expect_error(plot2Markers(.lm.obj = list()))
  expect_error(plotSamplePCA(.lm.obj = list()))
  expect_error(plotTradStats(.lm.obj = list(), .stats.obj = list()))
  expect_error(plotTradPerc(.lm.obj = list()))
  expect_error(plotAbundance(.lm.obj = list()))
  expect_error(plotDEA(.lm.obj = list(), .dea.obj = list(), .coefs = "test"))
})

test_that("interactFeatPlot validates input structure", {
  # Test that interactFeatPlot expects proper structure
  expect_error(interactFeatPlot(.lm.obj = list()))
  
  # Test with minimal structure that should work
  .lm.obj <- list(
    graph = list(uwot = list(embedding = matrix(runif(20), ncol=2))),
    pca = list(embed = matrix(runif(20), ncol=2)),
    lm = data.frame(A=runif(10), B=runif(10)),
    interact.plot = list(lm.features = list(html = rep("feature", 10)))
  )
  
  result <- interactFeatPlot(.lm.obj)
  
  # Check result type based on whether ggiraph is available
  if (requireNamespace("ggiraph", quietly = TRUE)) {
    expect_true(inherits(x = result, what = "girafe"))
  } else {
    expect_true(inherits(x = result, what = "ggplot"))
  }
})
