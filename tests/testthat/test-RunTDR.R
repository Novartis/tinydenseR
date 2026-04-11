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
# .get_sample_matrix tests
# ======================================================================

test_that(".get_sample_matrix returns identical matrix for files backend", {
  test_data <- create_test_lm_obj(n_cells = 10, n_markers = 3, n_samples = 2)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = c("marker_1", "marker_2", "marker_3"),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  result <- tinydenseR:::.get_sample_matrix(NULL, tdr, 1)
  expected <- readRDS(tdr@cells[[1]])
  expect_identical(result, expected)
})

# ======================================================================
# RunTDR.default tests
# ======================================================================

test_that("RunTDR.default runs full pipeline and returns TDRObj", {
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

  result <- RunTDR(tdr, .verbose = FALSE, .seed = 42)
  expect_true(is.TDRObj(result))
  expect_true(!is.null(result@density$fdens))
})

test_that("RunTDR.default with .celltype.vec runs celltyping", {
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

  # Build a named character vector of per-cell labels
  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  ct_vec <- stats::setNames(
    sample(c("TypeA", "TypeB"), length(all_cells), replace = TRUE),
    all_cells
  )

  result <- RunTDR(tdr, .celltype.vec = ct_vec, .verbose = FALSE, .seed = 42)
  expect_true(is.TDRObj(result))
  expect_true(!is.null(result@landmark.annot$celltyping$ids))
  # $mode no longer stored (removed in multi-solution redesign)
  expect_null(result@landmark.annot$celltyping$mode)
})

test_that("RunTDR.default errors on .celltype.vec conflict", {
  test_data <- create_test_lm_obj(n_cells = 10, n_markers = 3, n_samples = 2)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = c("marker_1", "marker_2", "marker_3"),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  # Set an existing celltype.vec in config
  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  old_vec <- stats::setNames(rep("Old", length(all_cells)), all_cells)
  tdr@config$celltype.vec <- old_vec

  new_vec <- stats::setNames(rep("New", length(all_cells)), all_cells)
  expect_error(
    RunTDR(tdr, .celltype.vec = new_vec, .verbose = FALSE),
    "overwrite"
  )
})

test_that("RunTDR.default with .celltype.vec.overwrite replaces vec", {
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

  # Set existing
  all_cells <- unlist(lapply(tdr@cells, function(p) rownames(readRDS(p))))
  old_vec <- stats::setNames(rep("Old", length(all_cells)), all_cells)
  tdr@config$celltype.vec <- old_vec

  new_vec <- stats::setNames(
    sample(c("TypeA", "TypeB"), length(all_cells), replace = TRUE),
    all_cells
  )
  # Should NOT error with .celltype.vec.overwrite = TRUE
  result <- RunTDR(tdr, .celltype.vec = new_vec,
                   .celltype.vec.overwrite = TRUE, .verbose = FALSE, .seed = 42)
  expect_true(is.TDRObj(result))
  expect_true(!is.null(result@landmark.annot$celltyping$ids))
})

# ======================================================================
# RunTDR.Seurat tests
# ======================================================================

test_that("RunTDR.Seurat stores TDRObj in Misc and returns Seurat", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  counts <- SingleCellExperiment::counts(sce)
  cell_meta <- data.frame(
    Sample = as.character(SummarizedExperiment::colData(sce)$Sample),
    row.names = colnames(sce)
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = counts, meta.data = cell_meta
  )

  result <- RunTDR(seurat_obj, .sample.var = "Sample",
                   .assay.type = "RNA", .verbose = FALSE, .seed = 42)

  expect_true(inherits(result, "Seurat"))
  tdr <- GetTDR(result)
  expect_true(is.TDRObj(tdr))
  expect_true(!is.null(tdr@density$fdens))
})

# ======================================================================
# RunTDR.SingleCellExperiment tests
# ======================================================================

test_that("RunTDR.SingleCellExperiment stores TDRObj and returns SCE", {
  skip_on_cran()
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  result <- RunTDR(sce, .sample.var = "Sample", .assay = "counts",
                   .assay.type = "RNA", .verbose = FALSE, .seed = 42)

  expect_true(inherits(result, "SingleCellExperiment"))
  tdr <- GetTDR(result)
  expect_true(is.TDRObj(tdr))
  expect_true(!is.null(tdr@density$fdens))
})

# ======================================================================
# Cross-backend equivalence tests
# ======================================================================

test_that("Seurat backend produces equivalent results to SCE", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  # --- SCE path ---
  sce_result <- RunTDR(sce, .sample.var = "Sample", .assay = "counts",
                       .assay.type = "RNA", .verbose = FALSE, .seed = 42)
  tdr_sce <- GetTDR(sce_result)

  # --- Seurat path ---
  counts <- SingleCellExperiment::counts(sce)
  cell_meta <- data.frame(
    Sample = as.character(SummarizedExperiment::colData(sce)$Sample),
    row.names = colnames(sce)
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = counts, meta.data = cell_meta
  )

  seurat_result <- RunTDR(seurat_obj, .sample.var = "Sample",
                          .assay.type = "RNA", .verbose = FALSE, .seed = 42)
  tdr_seurat <- GetTDR(seurat_result)

  # Compare density matrices
  expect_equal(as.matrix(tdr_sce@density$fdens),
               as.matrix(tdr_seurat@density$fdens),
               tolerance = 1e-10)
})

# ======================================================================
# GetTDR tests
# ======================================================================

test_that("GetTDR.default returns TDRObj as-is", {
  tdr <- TDRObj(config = list(assay.type = "RNA"))
  expect_identical(GetTDR(tdr), tdr)
})

test_that("GetTDR errors on non-TDR object", {
  expect_error(GetTDR(data.frame()), "Cannot extract")
})

test_that("GetTDR.Seurat errors when no TDRObj stored", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
  expect_error(GetTDR(seurat_obj), "No TDRObj found")
})

test_that("GetTDR.SingleCellExperiment errors when no TDRObj stored", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Matrix::sparseMatrix(
      i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
      dims = c(3, 3),
      dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
    ))
  )
  expect_error(GetTDR(sce), "No TDRObj found")
})

# ======================================================================
# RunTDR.HDF5AnnData and RunTDR.character tests
# ======================================================================

# Helper: create h5ad via Python anndata (anndataR-written h5ad files
# are not compatible with BPCells::open_matrix_anndata_hdf5).
.make_test_h5ad_runtdr <- function(file_path,
                                   n_cells = 200,
                                   n_genes = 50,
                                   n_samples = 4,
                                   use_layers_counts = FALSE) {

  python_script <- sprintf('
import anndata
import numpy as np
import scipy.sparse as sp
import pandas as pd

np.random.seed(42)
n_obs, n_var, n_samples = %d, %d, %d

counts = sp.random(n_obs, n_var, density=0.3, format="csr", dtype=np.float32)
counts.data = np.round(counts.data * 10).astype(np.float32)

cells_per_sample = [n_obs // n_samples] * n_samples
for i in range(n_obs %% n_samples):
    cells_per_sample[i] += 1

sample_ids = []
for i, c in enumerate(cells_per_sample):
    sample_ids.extend([f"sample{i+1}"] * c)

cell_ids = []
cell_counters = [0] * n_samples
for sid in sample_ids:
    idx = int(sid.replace("sample", "")) - 1
    cell_counters[idx] += 1
    cell_ids.append(f"{sid}_cell_{cell_counters[idx]}")

gene_ids = [f"gene_{i+1}" for i in range(n_var)]

obs = pd.DataFrame({
    "sample_id": pd.Categorical(sample_ids),
    "group": pd.Categorical(["A", "B"] * (n_obs // 2) + ["A"] * (n_obs %% 2))
}, index=cell_ids)

var = pd.DataFrame({
    "gene_name": gene_ids
}, index=gene_ids)

use_layers = %s

if use_layers:
    adata = anndata.AnnData(obs=obs, var=var)
    adata.layers["counts"] = counts.tocsc()
else:
    adata = anndata.AnnData(X=counts.tocsc(), obs=obs, var=var)

adata.write_h5ad("%s")
', n_cells, n_genes, n_samples,
  ifelse(use_layers_counts, "True", "False"),
  file_path)

  py_file <- tempfile(fileext = ".py")
  on.exit(unlink(py_file), add = TRUE)
  writeLines(python_script, py_file)

  result <- system2("python3", py_file, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop("Failed to create test h5ad fixture:\n",
         paste(result, collapse = "\n"))
  }

  invisible(file_path)
}

.has_python_anndata_runtdr <- function() {
  result <- tryCatch(
    system2("python3", c("-c", shQuote("import anndata")),
            stdout = TRUE, stderr = TRUE),
    error = function(e) "ERROR"
  )
  status <- attr(result, "status")
  is.null(status) || status == 0
}

test_that("RunTDR dispatches to HDF5AnnData method", {
  skip_if_not_installed("anndataR")
  # S3 dispatch reaches RunTDR.HDF5AnnData — verify it is registered
  expect_true(is.function(getS3method("RunTDR", "HDF5AnnData")))
})

test_that("RunTDR dispatches to character method", {
  expect_true(is.function(getS3method("RunTDR", "character")))
})

test_that("RunTDR.character runs full pipeline from h5ad file path", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata_runtdr(), "python3 with anndata not available")

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad_runtdr(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                         use_layers_counts = FALSE)

  result <- RunTDR(tmp_h5ad,
                   .sample.var = "sample_id",
                   .assay.type = "RNA",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE,
                   .seed = 42)

  expect_true(is.TDRObj(result))
  expect_true(!is.null(result@density$fdens))
  expect_equal(result@config$backend, "matrix")
  expect_equal(nrow(result@metadata), 4L)
})

test_that("RunTDR.character errors on missing file", {
  expect_error(RunTDR("nonexistent.h5ad", .sample.var = "x"),
               "File not found")
})

test_that("RunTDR.character errors on non-h5ad extension", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  writeLines("a,b", tmp)
  expect_error(RunTDR(tmp, .sample.var = "x"),
               "does not have an \\.h5ad extension")
})

test_that("RunTDR.character errors on invalid .sample.var", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata_runtdr(), "python3 with anndata not available")

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad_runtdr(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  expect_error(
    RunTDR(tmp_h5ad, .sample.var = "NonExistentCol", .verbose = FALSE),
    "not found in h5ad obs"
  )
})

test_that("RunTDR.HDF5AnnData runs full pipeline via S3 dispatch", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata_runtdr(), "python3 with anndata not available")

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad_runtdr(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                         use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  # Call via generic — S3 dispatch should find RunTDR.HDF5AnnData
  result <- RunTDR(x = x,
                   .sample.var = "sample_id",
                   .assay.type = "RNA",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE,
                   .seed = 42)

  expect_true(is.TDRObj(result))
  expect_true(!is.null(result@density$fdens))
  expect_equal(result@config$backend, "matrix")
})

test_that("RunTDR.HDF5AnnData errors on invalid .sample.var", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata_runtdr(), "python3 with anndata not available")

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")

  .make_test_h5ad_runtdr(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  expect_error(
    RunTDR(x = x,
           .sample.var = "NonExistentCol",
           .assay.type = "RNA",
           .verbose = FALSE),
    "not found in x\\$obs"
  )
})
