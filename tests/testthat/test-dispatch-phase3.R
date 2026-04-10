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
# Plot dispatch tests — verify S3 dispatch on TDRObj, Seurat, SCE
# ======================================================================

test_that("plotUMAP dispatches on TDRObj", {
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)
  p <- plotUMAP(tdr)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotPCA dispatches on TDRObj", {
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)
  p <- plotPCA(tdr)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotHeatmap dispatches on TDRObj", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("gridExtra")

  n_points <- 20
  n_clusters <- 3
  n_markers <- 4
  
  # Need full assay data for on-the-fly pheatmap computation
  expr_mat <- matrix(runif(n_points * n_markers), nrow = n_points, ncol = n_markers,
                     dimnames = list(paste0("lm", seq_len(n_points)),
                                    paste0("marker", seq_len(n_markers))))
  tdr <- TDRObj(
    config = list(assay.type = "cyto", markers = colnames(expr_mat)),
    landmark.embed = list(
      pca  = list(coord = matrix(
        runif(n_points * 2), ncol = 2,
        dimnames = list(NULL, c("PC1", "PC2"))
      )),
      umap = list(coord = matrix(runif(n_points * 2), ncol = 2)),
      le   = list()
    ),
    landmark.annot = list(
      clustering = list(
        ids = factor(sample(seq_len(n_clusters), n_points, replace = TRUE))
      )
    ),
    assay = list(expr = expr_mat),
    graphs = list(adj.matrix = NULL, snn = NULL, fgraph = NULL),
    density = list(),
    results = list()
  )

  p <- plotHeatmap(tdr)
  expect_true(inherits(p, "gtable") || inherits(p, "grob") ||
              inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotUMAP dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)
  p <- plotUMAP(seurat_obj)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotPCA dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)
  p <- plotPCA(seurat_obj)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotUMAP dispatches on SingleCellExperiment", {
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  sce <- SetTDR(sce, tdr)
  p <- plotUMAP(sce)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

# ======================================================================
# SetTDR round-trip for all container types
# ======================================================================

test_that("SetTDR + GetTDR round-trip on Seurat preserves TDRObj", {
  skip_if_not_installed("SeuratObject")
  tdr <- TDRObj(config = list(assay.type = "cyto"))

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)
  retrieved <- GetTDR(seurat_obj)
  expect_true(is.TDRObj(retrieved))
  expect_identical(retrieved@config$assay.type, "cyto")
})

test_that("SetTDR + GetTDR round-trip on SCE preserves TDRObj", {
  tdr <- TDRObj(config = list(assay.type = "cyto"))

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  sce <- SetTDR(sce, tdr)
  retrieved <- GetTDR(sce)
  expect_true(is.TDRObj(retrieved))
  expect_identical(retrieved@config$assay.type, "cyto")
})

# ======================================================================
# Round-trip regression: stepwise dispatch vs monolithic RunTDR
# ======================================================================

test_that("Stepwise dispatch matches RunTDR pipeline on TDRObj", {
  skip_on_cran()
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = c("marker_1", "marker_2", "marker_3", "marker_4", "marker_5"),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  # Monolithic
  result_mono <- RunTDR(tdr, .verbose = FALSE, .seed = 42)

  # Stepwise
  result_step <- tdr
  result_step <- get.landmarks(result_step, .verbose = FALSE, .seed = 42)
  result_step <- get.graph(result_step, .verbose = FALSE, .seed = 42)
  result_step <- get.map(result_step, .verbose = FALSE, .seed = 42)

  # Densities should match
  expect_equal(
    as.matrix(result_mono@density$fdens),
    as.matrix(result_step@density$fdens),
    tolerance = 1e-10
  )
})

# ======================================================================
# Generic existence tests
# ======================================================================

test_that("all plot generics exist and are functions", {
  plot_fns <- c(
    "plotPCA", "plotUMAP", "plotBeeswarm", "plot2Markers",
    "plotSamplePCA", "plotSampleEmbedding", "plotTradStats",
    "plotTradPerc", "plotDensity", "plotPbDE", "plotDEA",
    "plotMarkerDE", "plotHeatmap",
    "plotPlsD", "plotPlsDHeatmap"
  )
  for (fn_name in plot_fns) {
    fn <- get(fn_name, envir = asNamespace("tinydenseR"))
    expect_true(is.function(fn), info = paste(fn_name, "is not a function"))
  }
})

test_that("all plot .TDRObj methods exist", {
  plot_methods <- c(
    "plotPCA.TDRObj", "plotUMAP.TDRObj", "plotBeeswarm.TDRObj",
    "plot2Markers.TDRObj", "plotSamplePCA.TDRObj",
    "plotSampleEmbedding.TDRObj", "plotTradStats.TDRObj",
    "plotTradPerc.TDRObj", "plotDensity.TDRObj", "plotPbDE.TDRObj",
    "plotDEA.TDRObj", "plotMarkerDE.TDRObj", "plotHeatmap.TDRObj",
    "plotHeatmap.TDRObj",
    "plotPlsD.TDRObj", "plotPlsDHeatmap.TDRObj"
  )
  for (fn_name in plot_methods) {
    fn <- get(fn_name, envir = asNamespace("tinydenseR"))
    expect_true(is.function(fn), info = paste(fn_name, "is not a function"))
  }
})
