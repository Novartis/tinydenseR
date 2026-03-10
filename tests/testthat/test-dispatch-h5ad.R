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
# Helper: create a minimal mock HDF5AnnData-like R6 object for dispatch
# ======================================================================

# We cannot rely on anndataR being installed in CI, so we create a
# lightweight mock R6 object whose class vector includes "HDF5AnnData".
# This is sufficient for testing S3 dispatch through the thin wrappers.

create_mock_h5ad <- function(tdr = NULL) {
  obj <- list(uns = list())
  if (!is.null(tdr)) obj$uns[["tdr.obj"]] <- tdr
  class(obj) <- "HDF5AnnData"
  obj
}

# Provide `$` and `$<-` for our mock so uns assignment works like an R6 object
`$.HDF5AnnData` <- function(x, name) x[[name]]
`$<-.HDF5AnnData` <- function(x, name, value) { x[[name]] <- value; x }

# ======================================================================
# SetTDR / GetTDR round-trip for HDF5AnnData
# ======================================================================

test_that("SetTDR + GetTDR round-trip on HDF5AnnData preserves TDRObj", {
  tdr <- TDRObj(config = list(assay.type = "RNA"))
  h5ad <- create_mock_h5ad()

  h5ad <- SetTDR(h5ad, tdr)
  retrieved <- GetTDR(h5ad)
  expect_true(is.TDRObj(retrieved))
  expect_identical(retrieved@config$assay.type, "RNA")
})

test_that("GetTDR.HDF5AnnData errors when no TDRObj is stored", {
  h5ad <- create_mock_h5ad()
  expect_error(GetTDR(h5ad), "No TDRObj found")
})

# ======================================================================
# Pattern A dispatch tests (mutating, no .source)
# ======================================================================

test_that("get.graph dispatches on HDF5AnnData and returns HDF5AnnData", {
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

  h5ad <- create_mock_h5ad(tdr)
  result <- get.graph(h5ad, .k = 3, .verbose = FALSE, .seed = 42)
  expect_true(inherits(result, "HDF5AnnData"))
  tdr_out <- GetTDR(result)
  expect_true(!is.null(tdr_out@graphs$adj.matrix))
})

test_that("celltyping dispatches on HDF5AnnData and returns HDF5AnnData", {
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

  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  ct_vec <- stats::setNames(
    sample(c("TypeA", "TypeB"), length(all_cells), replace = TRUE),
    all_cells
  )

  h5ad <- create_mock_h5ad(tdr)
  result <- celltyping(h5ad, .celltyping.map = ct_vec, .verbose = FALSE)
  expect_true(inherits(result, "HDF5AnnData"))
  tdr_out <- GetTDR(result)
  expect_true(!is.null(tdr_out@landmark.annot$celltyping$ids))
})

# ======================================================================
# Pattern C dispatch tests (plot functions)
# ======================================================================

test_that("plotUMAP dispatches on HDF5AnnData", {
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)
  h5ad <- create_mock_h5ad(tdr)
  p <- plotUMAP(h5ad)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plotPCA dispatches on HDF5AnnData", {
  tdr <- create_mock_graph_obj(n_points = 20, n_clusters = 3)
  h5ad <- create_mock_h5ad(tdr)
  p <- plotPCA(h5ad)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

# ======================================================================
# Generic method existence tests
# ======================================================================

test_that("all HDF5AnnData dispatch methods exist", {
  stat_methods <- c(
    "get.lm.HDF5AnnData", "get.graph.HDF5AnnData",
    "celltyping.HDF5AnnData", "lm.cluster.HDF5AnnData",
    "get.features.HDF5AnnData", "get.embedding.HDF5AnnData",
    "get.specDE.HDF5AnnData", "get.nmfDE.HDF5AnnData",
    "get.plsDE.HDF5AnnData",
    "get.landmarks.HDF5AnnData", "get.map.HDF5AnnData",
    "get.pbDE.HDF5AnnData", "get.markerDE.HDF5AnnData",
    "goi.summary.HDF5AnnData"
  )
  for (fn_name in stat_methods) {
    fn <- get(fn_name, envir = asNamespace("tinydenseR"))
    expect_true(is.function(fn), info = paste(fn_name, "is not a function"))
  }
})

test_that("all HDF5AnnData plot dispatch methods exist", {
  plot_methods <- c(
    "plotPCA.HDF5AnnData", "plotUMAP.HDF5AnnData",
    "plotBeeswarm.HDF5AnnData", "plot2Markers.HDF5AnnData",
    "plotSamplePCA.HDF5AnnData", "plotSampleEmbedding.HDF5AnnData",
    "plotTradStats.HDF5AnnData", "plotTradPerc.HDF5AnnData",
    "plotDensity.HDF5AnnData", "plotPbDE.HDF5AnnData",
    "plotDEA.HDF5AnnData", "plotMarkerDE.HDF5AnnData",
    "plotHeatmap.HDF5AnnData", "plotSpecDE.HDF5AnnData",
    "plotSpecDEHeatmap.HDF5AnnData", "plotNmfDE.HDF5AnnData",
    "plotNmfDEHeatmap.HDF5AnnData", "plotPlsDE.HDF5AnnData",
    "plotPlsDEHeatmap.HDF5AnnData"
  )
  for (fn_name in plot_methods) {
    fn <- get(fn_name, envir = asNamespace("tinydenseR"))
    expect_true(is.function(fn), info = paste(fn_name, "is not a function"))
  }
})

test_that("GetTDR.HDF5AnnData and SetTDR.HDF5AnnData exist", {
  expect_true(is.function(
    get("GetTDR.HDF5AnnData", envir = asNamespace("tinydenseR"))
  ))
  expect_true(is.function(
    get("SetTDR.HDF5AnnData", envir = asNamespace("tinydenseR"))
  ))
})
