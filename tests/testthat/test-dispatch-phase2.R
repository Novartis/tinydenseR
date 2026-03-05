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
# Argument routing tests
# ======================================================================

test_that("RunTDR does not forward graph args to get.landmarks", {
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

  # .cl.resolution.parameter is a get.graph arg
  # This should NOT error — it should be routed only to get.graph
  expect_no_error(
    RunTDR(tdr, .cl.resolution.parameter = 2, .verbose = FALSE, .seed = 42)
  )
})

test_that("RunTDR warns on unknown arguments", {
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

  expect_warning(
    RunTDR(tdr, .bogus.param = 42, .verbose = FALSE, .seed = 42),
    "Unknown arguments"
  )
})

# ======================================================================
# Tier 2 dispatch tests
# ======================================================================

test_that("get.landmarks.TDRObj returns TDRObj with landmarks", {
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

  result <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)
  expect_true(is.TDRObj(result))
  expect_false(is.null(result@assay$raw))
})

# ======================================================================
# get.lm dispatch on Seurat
# ======================================================================

test_that("get.lm dispatches on Seurat and returns Seurat with lm results", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  # Use well-separated blobs so Leiden finds multiple clusters
  set.seed(42)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  blob2 <- matrix(rnorm(n * m, mean = 5, sd = 0.5), nrow = n, ncol = m)

  n_samples <- 4L
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

  tdr <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("marker_", 1:m),
    .assay.type = "cyto",
    .verbose = FALSE
  )
  tdr <- get.landmarks(tdr, .verbose = FALSE, .seed = 42)
  tdr <- get.graph(tdr, .k = 5, .scale = FALSE, .verbose = FALSE, .seed = 42)
  tdr <- get.map(tdr, .verbose = FALSE, .seed = 42)

  # Create design matrix: need as many rows as samples
  design <- cbind(intercept = 1,
                  group = as.numeric(tdr@metadata$group == "A"))

  # Store TDR in a Seurat

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_obj <- SetTDR(seurat_obj, tdr)

  result <- get.lm(seurat_obj, .design = design, .verbose = FALSE, .seed = 42)
  expect_true(inherits(result, "Seurat"))
  tdr_out <- GetTDR(result)
  expect_false(is.null(tdr_out@results$lm))
})
