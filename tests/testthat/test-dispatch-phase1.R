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

# ======================================================================
# SetTDR round-trip tests
# ======================================================================

test_that("SetTDR.default on a TDRObj returns the new TDRObj directly", {
  old_tdr <- TDRObj(config = list(assay.type = "RNA"))
  new_tdr <- TDRObj(config = list(assay.type = "cyto"))
  result <- SetTDR(old_tdr, new_tdr)
  expect_identical(result, new_tdr)
})

test_that("SetTDR.Seurat stores and GetTDR.Seurat retrieves identically", {
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  tdr <- TDRObj(config = list(assay.type = "RNA"))

  seurat_obj <- SetTDR(seurat_obj, tdr)
  retrieved <- GetTDR(seurat_obj)
  expect_true(is.TDRObj(retrieved))
  expect_identical(retrieved@config$assay.type, "RNA")
})

test_that("SetTDR.SingleCellExperiment stores and GetTDR.SingleCellExperiment retrieves identically", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  tdr <- TDRObj(config = list(assay.type = "RNA"))

  sce <- SetTDR(sce, tdr)
  retrieved <- GetTDR(sce)
  expect_true(is.TDRObj(retrieved))
  expect_identical(retrieved@config$assay.type, "RNA")
})

test_that("SetTDR.default errors on a data.frame", {
  expect_error(SetTDR(data.frame(), TDRObj()), "Cannot store")
})

# ======================================================================
# get.graph dispatch tests
# ======================================================================

test_that("get.graph on TDRObj returns TDRObj with adj.matrix populated", {
  skip_on_cran()
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)
  tdr <- get.graph(tdr, .k = 3, .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(tdr))
  expect_true(!is.null(tdr@graphs$adj.matrix))
})

test_that("get.graph dispatches on Seurat and returns Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)

  # Create minimal Seurat and store the TDRObj
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)

  result <- get.graph(seurat_obj, .k = 3, .verbose = FALSE, .seed = 42)
  expect_true(inherits(result, "Seurat"))
  tdr_out <- GetTDR(result)
  expect_true(!is.null(tdr_out@graphs$adj.matrix))
})

test_that("get.graph dispatches on SCE and returns SCE", {
  skip_on_cran()
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)

  # Create minimal SCE and store the TDRObj
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  sce <- SetTDR(sce, tdr)

  result <- get.graph(sce, .k = 3, .verbose = FALSE, .seed = 42)
  expect_true(inherits(result, "SingleCellExperiment"))
  tdr_out <- GetTDR(result)
  expect_true(!is.null(tdr_out@graphs$adj.matrix))
})

# ======================================================================
# celltyping dispatch test
# ======================================================================

test_that("celltyping dispatches on Seurat and returns Seurat with celltyping populated", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)
  tdr <- get.graph(tdr, .k = 3, .verbose = FALSE, .seed = 42)

  # Build a per-cell label vector matching ALL cells (not just landmarks)
  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  ct_vec <- stats::setNames(
    sample(c("TypeA", "TypeB"), length(all_cells), replace = TRUE),
    all_cells
  )

  # Store in Seurat
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)

  result <- celltyping(seurat_obj, .celltyping.map = ct_vec, .verbose = FALSE)
  expect_true(inherits(result, "Seurat"))
  tdr_out <- GetTDR(result)
  expect_true(!is.null(tdr_out@landmark.annot$celltyping$ids))
})
