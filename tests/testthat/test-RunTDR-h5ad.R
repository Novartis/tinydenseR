#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ======================================================================
# Test fixture helper
# ======================================================================

#' Create a minimal h5ad test fixture using Python anndata
#'
#' Writes a small h5ad with known values for round-trip testing.
#' Uses Python anndata via system call because anndataR-written h5ad
#' files are not compatible with BPCells::open_matrix_anndata_hdf5
#' (BPCells 0.3.x). Python-generated h5ad files work correctly.
#'
#' @param file_path Where to write the h5ad
#' @param n_cells Total number of cells
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param use_layers_counts If TRUE, store matrix in /layers/counts.
#'   If FALSE, store only in /X.
#' @return The file path (invisibly)
#' @noRd
.make_test_h5ad <- function(file_path,
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

# Helper to check if python3 with anndata is available
.has_python_anndata <- function() {
  result <- tryCatch(
    system2("python3", c("-c", shQuote("import anndata")),
            stdout = TRUE, stderr = TRUE),
    error = function(e) "ERROR"
  )
  status <- attr(result, "status")
  is.null(status) || status == 0
}

# ======================================================================
# .h5ad_resolve_counts_group tests
# ======================================================================

test_that(".h5ad_resolve_counts_group auto-detects /X", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2,
                  use_layers_counts = FALSE)

  result <- tinydenseR:::.h5ad_resolve_counts_group(tmp_h5ad, NULL)
  expect_equal(result, "/X")
})

test_that(".h5ad_resolve_counts_group auto-detects /layers/counts", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2,
                  use_layers_counts = TRUE)

  result <- tinydenseR:::.h5ad_resolve_counts_group(tmp_h5ad, NULL)
  expect_equal(result, "/layers/counts")
})

test_that(".h5ad_resolve_counts_group user override works", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2,
                  use_layers_counts = FALSE)

  result <- tinydenseR:::.h5ad_resolve_counts_group(tmp_h5ad, "/X")
  expect_equal(result, "/X")
})

test_that(".h5ad_resolve_counts_group errors on invalid group", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2,
                  use_layers_counts = FALSE)

  expect_error(
    tinydenseR:::.h5ad_resolve_counts_group(tmp_h5ad, "/nonexistent"),
    "could not be opened"
  )
})

test_that(".h5ad_resolve_counts_group errors when no valid group found", {
  skip_if_not_installed("BPCells")
  skip_on_cran()

  # Create a minimal HDF5 file with no valid matrix groups
  tmp_file <- withr::local_tempfile(fileext = ".h5ad")
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    rhdf5::h5createFile(tmp_file)
    rhdf5::h5createGroup(tmp_file, "obs")
    rhdf5::h5createGroup(tmp_file, "var")
    rhdf5::H5close()
  } else {
    skip("rhdf5 not available for creating minimal test fixture")
  }

  expect_error(
    tinydenseR:::.h5ad_resolve_counts_group(tmp_file, NULL),
    "Could not find a valid count matrix"
  )
})

# ======================================================================
# RunTDR.HDF5AnnData full pipeline tests
# ======================================================================

test_that("RunTDR.HDF5AnnData full pipeline with /X", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  expect_equal(result@config$backend, "matrix")
  expect_equal(nrow(result@metadata), 4L)
  expect_false(is.null(result@density$fdens))
  expect_true(all(paste0("sample", 1:4) %in% rownames(result@metadata)))
})

test_that("RunTDR.HDF5AnnData full pipeline with /layers/counts", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = TRUE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .h5ad.group = "/layers/counts",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  expect_equal(result@config$backend, "matrix")
  expect_equal(nrow(result@metadata), 4L)
  expect_false(is.null(result@density$fdens))
})

test_that("RunTDR.HDF5AnnData BPCells cache hit", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- file.path(withr::local_tempdir(), "bp_cache")

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  # First call: creates the BPCells directory
  result1 <- RunTDR(x,
                    .sample.var = "sample_id",
                    .bpcells.dir = bp_dir,
                    .nPC = 3,
                    .verbose = FALSE)

  expect_true(dir.exists(bp_dir))
  expect_s4_class(result1, "TDRObj")

  # Second call: should hit the cache
  x2 <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")
  output <- capture.output(
    result2 <- RunTDR(x2,
                      .sample.var = "sample_id",
                      .bpcells.dir = bp_dir,
                      .nPC = 3,
                      .verbose = TRUE),
    type = "output"
  )

  expect_true(any(grepl("cache hit", output, ignore.case = TRUE)))
  expect_s4_class(result2, "TDRObj")
})

test_that("RunTDR.HDF5AnnData with .bpcells.dir = NULL uses tempdir", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .bpcells.dir = NULL,
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
})

test_that("RunTDR.HDF5AnnData user-specified .h5ad.group override", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .h5ad.group = "/X",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
})

# ======================================================================
# Input validation tests
# ======================================================================

test_that("RunTDR.HDF5AnnData errors on missing .sample.var", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  expect_error(
    RunTDR(x, .sample.var = "nonexistent_column", .verbose = FALSE),
    "not found in x\\$obs"
  )
})

test_that("RunTDR.HDF5AnnData errors on non-string .sample.var", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  expect_error(
    RunTDR(x, .sample.var = 42, .verbose = FALSE),
    "must be a single character string"
  )
})

test_that("RunTDR.HDF5AnnData .min.cells.per.sample filtering", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  # Create fixture with uneven sample sizes:
  # sample1: 100 cells, sample2: 100 cells, sample3: 5 cells
  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  python_script <- sprintf('
import anndata
import numpy as np
import scipy.sparse as sp
import pandas as pd

np.random.seed(42)
n_var = 50
sample_sizes = [100, 100, 5]
n_obs = sum(sample_sizes)

counts = sp.random(n_obs, n_var, density=0.5, format="csr", dtype=np.float32)
counts.data = np.round(counts.data * 10 + 1).astype(np.float32)

sample_ids = []
cell_ids = []
for i, s in enumerate(sample_sizes):
    sid = f"sample{i+1}"
    for j in range(s):
        sample_ids.append(sid)
        cell_ids.append(f"{sid}_cell_{j+1}")

gene_ids = [f"gene_{i+1}" for i in range(n_var)]

obs = pd.DataFrame({
    "sample_id": pd.Categorical(sample_ids),
    "group": pd.Categorical(["A", "B"] * (n_obs // 2) + ["A"] * (n_obs %% 2))
}, index=cell_ids)
var = pd.DataFrame({"gene_name": gene_ids}, index=gene_ids)

adata = anndata.AnnData(X=counts.tocsc(), obs=obs, var=var)
adata.write_h5ad("%s")
', tmp_h5ad)

  py_file <- tempfile(fileext = ".py")
  on.exit(unlink(py_file), add = TRUE)
  writeLines(python_script, py_file)
  system2("python3", py_file)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .min.cells.per.sample = 10,
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  # sample3 has only 5 cells, so it should be excluded
  expect_equal(nrow(result@metadata), 2L)
  expect_false("sample3" %in% rownames(result@metadata))
})

# ======================================================================
# Data integrity tests
# ======================================================================

test_that("RunTDR.HDF5AnnData correct gene and cell names", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")
  expected_genes <- rownames(x$var)
  expected_cells <- rownames(x$obs)

  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  # Check gene names are in the source matrix
  mat <- result@config$source.env$mat
  expect_equal(rownames(mat), expected_genes)

  # Check that cell names in the result match expected cells
  result_cells <- unlist(lapply(result@cells, function(idx) colnames(mat)[idx]))
  expect_true(all(result_cells %in% expected_cells))
})

test_that("RunTDR.HDF5AnnData S3 dispatch works", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("anndataR")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  x <- anndataR::read_h5ad(tmp_h5ad, as = "HDF5AnnData")

  # Use RunTDR generic (not RunTDR.HDF5AnnData directly)
  result <- RunTDR(x,
                   .sample.var = "sample_id",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  expect_equal(result@config$backend, "matrix")
})

# ======================================================================
# Removal verification tests
# ======================================================================

test_that("h5ad backend case is removed from .get_sample_matrix", {
  # Check via the function body that no "h5ad" case exists in the switch
  fn_body <- deparse(tinydenseR:::.get_sample_matrix)
  h5ad_lines <- grep('"h5ad"', fn_body, value = TRUE)
  expect_length(h5ad_lines, 0L)
})

test_that("subset_h5ad is no longer exported", {
  exports <- getNamespaceExports("tinydenseR")
  expect_false("subset_h5ad" %in% exports)
})

# ======================================================================
# RunTDR.character tests (file-path-based h5ad support)
# ======================================================================

test_that("RunTDR.character full pipeline with /X", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  result <- RunTDR(tmp_h5ad,
                   .sample.var = "sample_id",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  expect_equal(result@config$backend, "matrix")
  expect_equal(nrow(result@metadata), 4L)
  expect_false(is.null(result@density$fdens))
  expect_true(all(paste0("sample", 1:4) %in% rownames(result@metadata)))
})

test_that("RunTDR.character full pipeline with /layers/counts", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = TRUE)

  result <- RunTDR(tmp_h5ad,
                   .sample.var = "sample_id",
                   .h5ad.group = "/layers/counts",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  expect_s4_class(result, "TDRObj")
  expect_equal(result@config$backend, "matrix")
  expect_equal(nrow(result@metadata), 4L)
})

test_that("RunTDR.character correct gene and cell names", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- withr::local_tempdir()

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  expected_genes <- tinydenseR:::.h5ad_read_var_names(tmp_h5ad)
  expected_cells <- tinydenseR:::.h5ad_read_obs_names(tmp_h5ad)

  result <- RunTDR(tmp_h5ad,
                   .sample.var = "sample_id",
                   .bpcells.dir = file.path(bp_dir, "bp"),
                   .nPC = 3,
                   .verbose = FALSE)

  mat <- result@config$source.env$mat
  expect_equal(rownames(mat), expected_genes)

  result_cells <- unlist(lapply(result@cells, function(idx) colnames(mat)[idx]))
  expect_true(all(result_cells %in% expected_cells))
})

test_that("RunTDR.character BPCells cache hit", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  bp_dir <- file.path(withr::local_tempdir(), "bp_cache")

  .make_test_h5ad(tmp_h5ad, n_cells = 200, n_genes = 50, n_samples = 4,
                  use_layers_counts = FALSE)

  # First call: creates the BPCells directory
  result1 <- RunTDR(tmp_h5ad,
                    .sample.var = "sample_id",
                    .bpcells.dir = bp_dir,
                    .nPC = 3,
                    .verbose = FALSE)

  expect_true(dir.exists(bp_dir))
  expect_s4_class(result1, "TDRObj")

  # Second call: should hit the cache
  output <- capture.output(
    result2 <- RunTDR(tmp_h5ad,
                      .sample.var = "sample_id",
                      .bpcells.dir = bp_dir,
                      .nPC = 3,
                      .verbose = TRUE),
    type = "output"
  )

  expect_true(any(grepl("cache hit", output, ignore.case = TRUE)))
  expect_s4_class(result2, "TDRObj")
})

test_that("RunTDR.character errors on missing .sample.var", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  expect_error(
    RunTDR(tmp_h5ad, .sample.var = "nonexistent_column", .verbose = FALSE),
    "not found in h5ad obs"
  )
})

test_that("RunTDR.character errors on non-string .sample.var", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("rhdf5")
  skip_if(!.has_python_anndata(), "python3 with anndata not available")
  skip_on_cran()

  tmp_h5ad <- withr::local_tempfile(fileext = ".h5ad")
  .make_test_h5ad(tmp_h5ad, n_cells = 40, n_genes = 10, n_samples = 2)

  expect_error(
    RunTDR(tmp_h5ad, .sample.var = 42, .verbose = FALSE),
    "must be a single character string"
  )
})

test_that("HDF5AnnData dispatch methods no longer exist", {
  # GetTDR and SetTDR .HDF5AnnData methods were removed
  expect_null(tryCatch(
    getS3method("GetTDR", "HDF5AnnData"),
    error = function(e) NULL
  ))
  expect_null(tryCatch(
    getS3method("SetTDR", "HDF5AnnData"),
    error = function(e) NULL
  ))
  # Spot-check: downstream dispatch methods also removed
  expect_null(tryCatch(
    getS3method("get.graph", "HDF5AnnData"),
    error = function(e) NULL
  ))
  expect_null(tryCatch(
    getS3method("plotUMAP", "HDF5AnnData"),
    error = function(e) NULL
  ))
})
