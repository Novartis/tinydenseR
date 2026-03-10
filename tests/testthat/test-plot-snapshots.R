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

# Helper: create a reproducible mock TDRObj with correct UMAP colnames
make_snapshot_obj <- function(n_points = 20, n_clusters = 3, seed = 42) {
  set.seed(seed)
  tdr <- create_mock_graph_obj(n_points = n_points, n_clusters = n_clusters)
  colnames(tdr@landmark.embed$umap$coord) <- c("umap.1", "umap.2")
  tdr
}

# ======================================================================
# plotPCA snapshots
# ======================================================================

test_that("plotPCA snapshot is consistent", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")

  tdr <- make_snapshot_obj()
  p <- plotPCA(tdr, .seed = 42)
  vdiffr::expect_doppelganger("plotPCA-base", p)
})

test_that("plotPCA with numeric feature snapshot is consistent", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")

  tdr <- make_snapshot_obj()
  numeric_feature <- runif(20) # same RNG state after make_snapshot_obj(seed=42)
  set.seed(42) # reset for deterministic feature
  numeric_feature <- runif(20)
  p <- plotPCA(tdr, .feature = numeric_feature, .color.label = "value", .seed = 42)
  vdiffr::expect_doppelganger("plotPCA-numeric-feature", p)
})

test_that("plotPCA with categorical feature snapshot is consistent", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")

  tdr <- make_snapshot_obj()
  cat_feature <- factor(rep(c("TypeA", "TypeB", "TypeC", "TypeD"), length.out = 20))
  p <- plotPCA(tdr, .feature = cat_feature, .color.label = "cell type", .seed = 42)
  vdiffr::expect_doppelganger("plotPCA-categorical-feature", p)
})

# ======================================================================
# plotUMAP snapshots
# ======================================================================

test_that("plotUMAP snapshot is consistent", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")

  tdr <- make_snapshot_obj()
  p <- plotUMAP(tdr, .seed = 42)
  vdiffr::expect_doppelganger("plotUMAP-base", p)
})

# ======================================================================
# Cross-backend equivalence: Seurat
# ======================================================================

test_that("plotPCA produces same plot data from Seurat wrapper", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("SeuratObject")

  tdr <- make_snapshot_obj()

  p_direct <- plotPCA(tdr, .seed = 42)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)
  p_seurat <- plotPCA(seurat_obj, .seed = 42)

  expect_identical(p_direct$data, p_seurat$data)
  vdiffr::expect_doppelganger("plotPCA-seurat", p_seurat)
})

test_that("plotUMAP produces same plot data from Seurat wrapper", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("SeuratObject")

  tdr <- make_snapshot_obj()

  p_direct <- plotUMAP(tdr, .seed = 42)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)
  p_seurat <- plotUMAP(seurat_obj, .seed = 42)

  expect_identical(p_direct$data, p_seurat$data)
  vdiffr::expect_doppelganger("plotUMAP-seurat", p_seurat)
})

# ======================================================================
# Cross-backend equivalence: SingleCellExperiment
# ======================================================================

test_that("plotPCA produces same plot data from SCE wrapper", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("SingleCellExperiment")

  tdr <- make_snapshot_obj()

  p_direct <- plotPCA(tdr, .seed = 42)

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  sce <- SetTDR(sce, tdr)
  p_sce <- plotPCA(sce, .seed = 42)

  expect_identical(p_direct$data, p_sce$data)
  vdiffr::expect_doppelganger("plotPCA-sce", p_sce)
})

test_that("plotUMAP produces same plot data from SCE wrapper", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("SingleCellExperiment")

  tdr <- make_snapshot_obj()

  p_direct <- plotUMAP(tdr, .seed = 42)

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  sce <- SetTDR(sce, tdr)
  p_sce <- plotUMAP(sce, .seed = 42)

  expect_identical(p_direct$data, p_sce$data)
  vdiffr::expect_doppelganger("plotUMAP-sce", p_sce)
})
