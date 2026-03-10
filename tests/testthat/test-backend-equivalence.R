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
# Helpers: wrap a TDRObj in Seurat or SCE containers
# ======================================================================

wrap_in_seurat <- function(tdr) {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  SetTDR(seurat_obj, tdr)
}

wrap_in_sce <- function(tdr) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  SetTDR(sce, tdr)
}

# ======================================================================
# get.lm equivalence: bare TDRObj vs Seurat vs SCE
# ======================================================================

test_that("get.lm produces identical results via TDRObj, Seurat, and SCE", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  # Use well-separated blobs so Leiden finds multiple clusters
  set.seed(42)
  n <- 200; m <- 5; n_samples <- 4L
  blob1 <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  blob2 <- matrix(rnorm(n * m, mean = 5, sd = 0.5), nrow = n, ncol = m)

  .cells <- lapply(seq_len(n_samples), function(i) {
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    dimnames(mat) <- list(
      paste0("s", i, "_", seq_len(nrow(mat))),
      paste0("marker_", seq_len(m))
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    return(uri)
  })
  names(.cells) <- paste0("sample", seq_len(n_samples))
  .meta <- data.frame(
    row.names = names(.cells),
    group = rep(c("A", "B"), length.out = n_samples)
  )
  on.exit(lapply(.cells, unlink), add = TRUE)

  # Build pipeline up to map stage
  tdr_base <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("marker_", 1:m),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr_base <- get.landmarks(tdr_base, .verbose = FALSE, .seed = 42)
  tdr_base <- get.graph(tdr_base, .k = 5, .scale = FALSE,
                        .verbose = FALSE, .seed = 42)
  tdr_base <- get.map(tdr_base, .verbose = FALSE, .seed = 42)
  design <- model.matrix(~ group, data = .meta)


  # --- bare TDRObj ---
  tdr_direct <- get.lm(tdr_base, .design = design, .verbose = FALSE, .seed = 42)
  res_direct <- tdr_direct@results$lm$default

  # --- Seurat wrapper ---
  seurat_obj <- wrap_in_seurat(tdr_base)
  seurat_result <- get.lm(seurat_obj, .design = design, .verbose = FALSE, .seed = 42)
  res_seurat <- GetTDR(seurat_result)@results$lm$default

  # --- SCE wrapper ---
  sce_obj <- wrap_in_sce(tdr_base)
  sce_result <- get.lm(sce_obj, .design = design, .verbose = FALSE, .seed = 42)
  res_sce <- GetTDR(sce_result)@results$lm$default

  # Coefficients

  expect_equal(res_direct$fit$coefficients, res_seurat$fit$coefficients,
               tolerance = 1e-12)
  expect_equal(res_direct$fit$coefficients, res_sce$fit$coefficients,
               tolerance = 1e-12)

  # P-values
  expect_equal(res_direct$fit$p.value, res_seurat$fit$p.value,
               tolerance = 1e-12)
  expect_equal(res_direct$fit$p.value, res_sce$fit$p.value,
               tolerance = 1e-12)

  # PCA-weighted q-values
  expect_equal(res_direct$fit$pca.weighted.q, res_seurat$fit$pca.weighted.q,
               tolerance = 1e-12)
  expect_equal(res_direct$fit$pca.weighted.q, res_sce$fit$pca.weighted.q,
               tolerance = 1e-12)

  # Traditional clustering coefficients
  expect_equal(res_direct$trad$clustering$fit$coefficients,
               res_seurat$trad$clustering$fit$coefficients,
               tolerance = 1e-12)
  expect_equal(res_direct$trad$clustering$fit$coefficients,
               res_sce$trad$clustering$fit$coefficients,
               tolerance = 1e-12)
})

# ======================================================================
# get.graph equivalence: bare TDRObj vs Seurat vs SCE
# ======================================================================

test_that("get.graph produces identical results via TDRObj and Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  # Build up to landmarks
  tdr_lm <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr_lm <- get.landmarks(tdr_lm, .verbose = FALSE, .seed = 42)

  # --- bare TDRObj ---
  tdr_graph <- get.graph(tdr_lm, .k = 3, .verbose = FALSE, .seed = 42)

  # --- Seurat wrapper ---
  seurat_obj <- wrap_in_seurat(tdr_lm)
  seurat_result <- get.graph(seurat_obj, .k = 3, .verbose = FALSE, .seed = 42)
  tdr_seurat <- GetTDR(seurat_result)

  # Adjacency matrix
  expect_identical(tdr_graph@graphs$adj.matrix, tdr_seurat@graphs$adj.matrix)

  # Clustering ids (same seed → identical clusters)
  expect_identical(tdr_graph@landmark.annot$clustering$ids,
                   tdr_seurat@landmark.annot$clustering$ids)
})

test_that("get.graph produces identical results via TDRObj and SCE", {
  skip_on_cran()

  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr_lm <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr_lm <- get.landmarks(tdr_lm, .verbose = FALSE, .seed = 42)

  # --- bare TDRObj ---
  tdr_graph <- get.graph(tdr_lm, .k = 3, .verbose = FALSE, .seed = 42)

  # --- SCE wrapper ---
  sce_obj <- wrap_in_sce(tdr_lm)
  sce_result <- get.graph(sce_obj, .k = 3, .verbose = FALSE, .seed = 42)
  tdr_sce <- GetTDR(sce_result)

  # Adjacency matrix
  expect_identical(tdr_graph@graphs$adj.matrix, tdr_sce@graphs$adj.matrix)

  # Clustering ids
  expect_identical(tdr_graph@landmark.annot$clustering$ids,
                   tdr_sce@landmark.annot$clustering$ids)
})

# ======================================================================
# get.features equivalence: bare TDRObj vs Seurat
# ======================================================================

test_that("get.features produces identical results via TDRObj and Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr_lm <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr_lm <- get.landmarks(tdr_lm, .verbose = FALSE, .seed = 42)

  # --- bare TDRObj ---
  tdr_feat <- get.features(tdr_lm)

  # --- Seurat wrapper ---
  seurat_obj <- wrap_in_seurat(tdr_lm)
  seurat_result <- get.features(seurat_obj)
  tdr_seurat <- GetTDR(seurat_result)

  # Feature results
  expect_identical(tdr_feat@results$features$lm.features$res,
                   tdr_seurat@results$features$lm.features$res)
})

# ======================================================================
# plotPCA data equivalence: bare TDRObj vs Seurat vs SCE
# ======================================================================

test_that("plotPCA produces identical plot data via TDRObj and Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  p_direct <- plotPCA(tdr)

  seurat_obj <- wrap_in_seurat(tdr)
  p_seurat <- plotPCA(seurat_obj)

  expect_identical(p_direct$data, p_seurat$data)
})

test_that("plotPCA produces identical plot data via TDRObj and SCE", {
  skip_on_cran()

  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  p_direct <- plotPCA(tdr)

  sce_obj <- wrap_in_sce(tdr)
  p_sce <- plotPCA(sce_obj)

  expect_identical(p_direct$data, p_sce$data)
})

test_that("plotUMAP produces identical plot data via TDRObj and Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  p_direct <- plotUMAP(tdr)

  seurat_obj <- wrap_in_seurat(tdr)
  p_seurat <- plotUMAP(seurat_obj)

  expect_identical(p_direct$data, p_seurat$data)
})

test_that("plotUMAP produces identical plot data via TDRObj and SCE", {
  skip_on_cran()

  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  p_direct <- plotUMAP(tdr)

  sce_obj <- wrap_in_sce(tdr)
  p_sce <- plotUMAP(sce_obj)

  expect_identical(p_direct$data, p_sce$data)
})
