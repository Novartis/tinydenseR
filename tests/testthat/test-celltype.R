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
  .tdr.obj <- TDRObj(
    graph = list(clustering = list(ids = factor(c("cluster.01", "cluster.02"))))
  )
  
  # Test cluster mapping to multiple celltypes
  .celltyping.map.multi <- list(
    "T_cells" = c("cluster.01"),
    "B_cells" = c("cluster.01")  # Same cluster mapped to multiple celltypes
  )
  expect_error(celltyping(.tdr.obj, .celltyping.map = .celltyping.map.multi),
               "Cluster\\(s\\) mapped to multiple cell types")
  
  # Test unmapped clusters
  .celltyping.map.incomplete <- list(
    "T_cells" = c("cluster.01")
    # cluster.02 is not mapped
  )
  expect_error(celltyping(.tdr.obj, .celltyping.map = .celltyping.map.incomplete),
               "Every cluster must be mapped to a cell type")
})

test_that("celltyping requires named list", {
  .tdr.obj <- TDRObj(
    graph = list(clustering = list(ids = factor(c("cluster.01"))))
  )
  
  # Test unnamed list
  .celltyping.map.unnamed <- list(c("cluster.01"))
  expect_error(celltyping(.tdr.obj, .celltyping.map = .celltyping.map.unnamed),
               "Cell type names are missing")
})

# ======================================================================
# Regression: categorical metadata without .celltype.vec (GH bug)
# ======================================================================

test_that("get.map succeeds when celltyping list exists but $ids is NULL", {
  skip_on_cran()

  # Build a TDRObj via the cyto pipeline (no Seurat dep needed)
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- run_full_pipeline(test_data, seed = 42, verbose = FALSE)

  # Simulate what import_cell_annotations does when .celltype.vec is NULL:
  # it stores named solutions but sets $ids back to NULL.
  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  fake_labels <- stats::setNames(
    sample(c("TypeA", "TypeB"), length(all_cells), replace = TRUE),
    all_cells
  )
  tdr <- celltyping(tdr, .celltyping.map = fake_labels,
                    .name = "Phase", .verbose = FALSE)

  # Confirm $ids is currently set

  expect_false(is.null(tdr@landmark.annot$celltyping$ids))

  # Now NULL out $ids (as import_cell_annotations does when no .celltype.vec)
  tdr@landmark.annot$celltyping$ids <- NULL

  # celltyping list is non-NULL but $ids is NULL — this used to crash get.map
  expect_false(is.null(tdr@landmark.annot$celltyping))
  expect_true(is.null(tdr@landmark.annot$celltyping$ids))

  # Also clear cellmap-level celltyping from the prior run
  tdr@cellmap$celltyping$ids <- NULL

  # Re-run get.map — must not error
  expect_no_error(
    tdr <- get.map(tdr, .verbose = FALSE, .seed = 42)
  )

  # Clustering should still work
  expect_false(is.null(tdr@cellmap$clustering$ids))

  # Stored solution should still be accessible
  expect_true("Phase" %in% list_celltyping_solutions(tdr))
})
