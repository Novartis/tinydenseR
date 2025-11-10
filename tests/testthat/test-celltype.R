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

# Test for celltyping

test_that("celltyping validates input mapping", {
  .lm.obj <- list(
    graph = list(clustering = list(ids = factor(c("cluster.01", "cluster.02"))))
  )
  
  # Test cluster mapping to multiple celltypes
  .celltyping.map.multi <- list(
    "T_cells" = c("cluster.01"),
    "B_cells" = c("cluster.01")  # Same cluster mapped to multiple celltypes
  )
  expect_error(celltyping(.lm.obj = .lm.obj, .celltyping.map = .celltyping.map.multi),
               "Cluster\\(s\\) mapped to multiple cell types")
  
  # Test unmapped clusters
  .celltyping.map.incomplete <- list(
    "T_cells" = c("cluster.01")
    # cluster.02 is not mapped
  )
  expect_error(celltyping(.lm.obj = .lm.obj, .celltyping.map = .celltyping.map.incomplete),
               "Every cluster must be mapped to a cell type")
})

test_that("celltyping requires named list", {
  .lm.obj <- list(
    graph = list(clustering = list(ids = factor(c("cluster.01"))))
  )
  
  # Test unnamed list
  .celltyping.map.unnamed <- list(c("cluster.01"))
  expect_error(celltyping(.lm.obj = .lm.obj, .celltyping.map = .celltyping.map.unnamed),
               "Cell type names are missing")
})
