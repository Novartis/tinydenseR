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

# ──────────────────────────────────────────────────────────────────────
# Test helpers: write CSR/CSC from a dense matrix into an HDF5 group
# ──────────────────────────────────────────────────────────────────────

#' Write a dense matrix as CSR into an existing HDF5 group
#' @param file HDF5 file path
#' @param group_path HDF5 group path (must already exist)
#' @param mat Dense matrix (n_obs x n_var)
#' @noRd
.write_csr <- function(file, group_path, mat) {
  n_obs <- nrow(mat)
  n_var <- ncol(mat)

  data    <- numeric(0)
  indices <- integer(0)
  indptr  <- integer(n_obs + 1L)
  indptr[1L] <- 0L

  for (i in seq_len(n_obs)) {
    nz <- which(mat[i, ] != 0)
    if (length(nz) > 0L) {
      data    <- c(data, as.numeric(mat[i, nz]))
      indices <- c(indices, as.integer(nz - 1L))  # 0-based
    }
    indptr[i + 1L] <- length(data)
  }

  indptr <- as.integer(indptr)

  rhdf5::h5write(data,    file, paste0(group_path, "/data"))
  rhdf5::h5write(indices, file, paste0(group_path, "/indices"))
  rhdf5::h5write(indptr,  file, paste0(group_path, "/indptr"))
  rhdf5::h5writeAttribute(c(n_obs, n_var), file, "shape",
                          h5loc = group_path)
  rhdf5::h5writeAttribute("csr_matrix", file, "encoding-type",
                          h5loc = group_path, asScalar = TRUE)
  rhdf5::h5writeAttribute("0.1.0", file, "encoding-version",
                          h5loc = group_path, asScalar = TRUE)
}

#' Write a dense matrix as CSC into an existing HDF5 group
#' @param file HDF5 file path
#' @param group_path HDF5 group path (must already exist)
#' @param mat Dense matrix (n_obs x n_var)
#' @noRd
.write_csc <- function(file, group_path, mat) {
  n_obs <- nrow(mat)
  n_var <- ncol(mat)

  data    <- numeric(0)
  indices <- integer(0)
  indptr  <- integer(n_var + 1L)
  indptr[1L] <- 0L

  for (j in seq_len(n_var)) {
    nz <- which(mat[, j] != 0)
    if (length(nz) > 0L) {
      data    <- c(data, as.numeric(mat[nz, j]))
      indices <- c(indices, as.integer(nz - 1L))  # 0-based
    }
    indptr[j + 1L] <- length(data)
  }

  indptr <- as.integer(indptr)

  rhdf5::h5write(data,    file, paste0(group_path, "/data"))
  rhdf5::h5write(indices, file, paste0(group_path, "/indices"))
  rhdf5::h5write(indptr,  file, paste0(group_path, "/indptr"))
  rhdf5::h5writeAttribute(c(n_obs, n_var), file, "shape",
                          h5loc = group_path)
  rhdf5::h5writeAttribute("csc_matrix", file, "encoding-type",
                          h5loc = group_path, asScalar = TRUE)
  rhdf5::h5writeAttribute("0.1.0", file, "encoding-version",
                          h5loc = group_path, asScalar = TRUE)
}

#' Create a minimal valid .h5ad file for testing
#' @param path File path to create
#' @param n_obs Number of cells
#' @param n_var Number of genes
#' @param encoding One of "csr_matrix", "csc_matrix", "array"
#' @param n_layers Number of additional layers (0 = none)
#' @param with_obsm Logical — add an X_pca embedding?
#' @return Invisible path
#' @noRd
create_test_h5ad <- function(path,
                             n_obs    = 100L,
                             n_var    = 50L,
                             encoding = "csr_matrix",
                             n_layers = 0L,
                             with_obsm = TRUE) {

  rhdf5::h5createFile(path)

  # Root attributes
  rhdf5::h5writeAttribute("anndata", path, "encoding-type",
                          h5loc = "/", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.1.0", path, "encoding-version",
                          h5loc = "/", asScalar = TRUE)

  # --- /X ---
  set.seed(42)
  vals <- rpois(n_obs * n_var, lambda = 0.3)
  dense_mat <- matrix(vals, nrow = n_obs, ncol = n_var)

  if (encoding == "array") {
    rhdf5::h5write(dense_mat, path, "/X")
    rhdf5::h5writeAttribute("array", path, "encoding-type",
                            h5loc = "/X", asScalar = TRUE)
    rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                            h5loc = "/X", asScalar = TRUE)
  } else if (encoding == "csr_matrix") {
    rhdf5::h5createGroup(path, "/X")
    .write_csr(path, "/X", dense_mat)
  } else if (encoding == "csc_matrix") {
    rhdf5::h5createGroup(path, "/X")
    .write_csc(path, "/X", dense_mat)
  }

  # --- /obs ---
  rhdf5::h5createGroup(path, "/obs")
  cell_names <- paste0("cell_", seq_len(n_obs))
  sample_ids <- paste0("sample_", rep(1:5, length.out = n_obs))

  rhdf5::h5write(cell_names, path, "/obs/_index")
  rhdf5::h5writeAttribute("array", path, "encoding-type",
                          h5loc = "/obs/_index", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/obs/_index", asScalar = TRUE)
  rhdf5::h5write(sample_ids, path, "/obs/sample_id")
  rhdf5::h5writeAttribute("array", path, "encoding-type",
                          h5loc = "/obs/sample_id", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/obs/sample_id", asScalar = TRUE)

  # Categorical column
  rhdf5::h5createGroup(path, "/obs/cell_type")
  categories <- c("T_cell", "B_cell", "Monocyte")
  set.seed(123)
  codes <- sample(0:2, n_obs, replace = TRUE)
  rhdf5::h5write(codes,      path, "/obs/cell_type/codes")
  rhdf5::h5write(categories, path, "/obs/cell_type/categories")
  rhdf5::h5writeAttribute("categorical", path, "encoding-type",
                          h5loc = "/obs/cell_type", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/obs/cell_type", asScalar = TRUE)
  rhdf5::h5writeAttribute(FALSE, path, "ordered",
                          h5loc = "/obs/cell_type")

  rhdf5::h5writeAttribute("dataframe", path, "encoding-type",
                          h5loc = "/obs", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/obs", asScalar = TRUE)
  rhdf5::h5writeAttribute("_index", path, "_index",
                          h5loc = "/obs", asScalar = TRUE)
  rhdf5::h5writeAttribute(c("sample_id", "cell_type"), path, "column-order",
                          h5loc = "/obs")

  # --- /var ---
  rhdf5::h5createGroup(path, "/var")
  gene_names <- paste0("gene_", seq_len(n_var))
  rhdf5::h5write(gene_names, path, "/var/_index")
  rhdf5::h5writeAttribute("array", path, "encoding-type",
                          h5loc = "/var/_index", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/var/_index", asScalar = TRUE)
  rhdf5::h5writeAttribute("dataframe", path, "encoding-type",
                          h5loc = "/var", asScalar = TRUE)
  rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                          h5loc = "/var", asScalar = TRUE)
  rhdf5::h5writeAttribute("_index", path, "_index",
                          h5loc = "/var", asScalar = TRUE)
  # rhdf5 cannot write character(0) as an attribute; use low-level API
  fid <- rhdf5::H5Fopen(path)
  gid <- rhdf5::H5Gopen(fid, "/var")
  tid <- rhdf5::H5Tcopy("H5T_C_S1")
  rhdf5::H5Tset_size(tid, 1L)
  sid <- rhdf5::H5Screate_simple(0)
  aid <- rhdf5::H5Acreate(gid, "column-order", tid, sid)
  rhdf5::H5Aclose(aid)
  rhdf5::H5Sclose(sid)
  rhdf5::H5Gclose(gid)
  rhdf5::H5Fclose(fid)

  # --- /obsm ---
  if (with_obsm) {
    rhdf5::h5createGroup(path, "/obsm")
    set.seed(99)
    pca_mat <- matrix(rnorm(n_obs * 10), nrow = n_obs, ncol = 10)
    rhdf5::h5write(pca_mat, path, "/obsm/X_pca")
    rhdf5::h5writeAttribute("array", path, "encoding-type",
                            h5loc = "/obsm/X_pca", asScalar = TRUE)
    rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                            h5loc = "/obsm/X_pca", asScalar = TRUE)
  }

  # --- /uns ---
  rhdf5::h5createGroup(path, "/uns")
  rhdf5::h5write("test_dataset", path, "/uns/dataset_name")

  # --- /layers ---
  if (n_layers > 0L) {
    rhdf5::h5createGroup(path, "/layers")
    for (li in seq_len(n_layers)) {
      layer_name <- paste0("layer_", li)
      set.seed(200 + li)
      layer_mat <- matrix(runif(n_obs * n_var), nrow = n_obs, ncol = n_var)
      layer_path <- paste0("/layers/", layer_name)
      rhdf5::h5write(layer_mat, path, layer_path)
      rhdf5::h5writeAttribute("array", path, "encoding-type",
                              h5loc = layer_path, asScalar = TRUE)
      rhdf5::h5writeAttribute("0.2.0", path, "encoding-version",
                              h5loc = layer_path, asScalar = TRUE)
    }
  }

  invisible(path)
}

#' Reconstruct a dense matrix from CSR components read from h5ad
#' @noRd
.csr_to_dense <- function(data, indices, indptr, n_row, n_col) {
  mat <- matrix(0, nrow = n_row, ncol = n_col)
  for (i in seq_len(n_row)) {
    start <- indptr[i] + 1L
    end   <- indptr[i + 1L]
    if (end >= start) {
      cols <- indices[start:end] + 1L
      mat[i, cols] <- data[start:end]
    }
  }
  mat
}

#' Reconstruct a dense matrix from CSC components read from h5ad
#' @noRd
.csc_to_dense <- function(data, indices, indptr, n_row, n_col) {
  mat <- matrix(0, nrow = n_row, ncol = n_col)
  for (j in seq_len(n_col)) {
    start <- indptr[j] + 1L
    end   <- indptr[j + 1L]
    if (end >= start) {
      rows <- indices[start:end] + 1L
      mat[rows, j] <- data[start:end]
    }
  }
  mat
}

#' Read the output matrix from a subset h5ad file as dense
#' @noRd
.read_output_matrix_as_dense <- function(path, n_row, n_col) {
  enc <- rhdf5::h5readAttributes(path, "/X")[["encoding-type"]]
  if (enc == "array") {
    return(rhdf5::h5read(path, "/X"))
  }
  data    <- rhdf5::h5read(path, "/X/data")
  indices <- rhdf5::h5read(path, "/X/indices")
  indptr  <- rhdf5::h5read(path, "/X/indptr")
  if (enc == "csr_matrix") {
    .csr_to_dense(data, indices, indptr, n_row, n_col)
  } else {
    .csc_to_dense(data, indices, indptr, n_row, n_col)
  }
}

# ──────────────────────────────────────────────────────────────────────
# Tests
# ──────────────────────────────────────────────────────────────────────

test_that("subset_h5ad works with CSR matrix", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  cell_idx <- c(10L, 5L, 50L, 75L, 1L)
  gene_idx <- c(3L, 1L, 20L, 40L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = gene_idx,
    .chunk.size  = 2L,
    .verbose     = FALSE
  )

  expect_true(file.exists(dst))

  # Verify dimensions
  obs_index <- rhdf5::h5read(dst, "/obs/_index")
  expect_equal(length(obs_index), length(cell_idx))

  var_index <- rhdf5::h5read(dst, "/var/_index")
  expect_equal(length(var_index), length(gene_idx))

  # Verify X values
  set.seed(42)
  original_dense <- matrix(rpois(100 * 50, lambda = 0.3), nrow = 100, ncol = 50)
  expected_sub <- original_dense[cell_idx, gene_idx]

  result_dense <- .read_output_matrix_as_dense(dst, length(cell_idx),
                                               length(gene_idx))
  expect_equal(result_dense, expected_sub)

  # Verify cell names are properly reordered
  src_cell_names <- paste0("cell_", seq_len(100))
  expect_equal(c(obs_index), src_cell_names[cell_idx])

  # Verify gene names are properly subset
  src_gene_names <- paste0("gene_", seq_len(50))
  expect_equal(c(var_index), src_gene_names[gene_idx])
})

test_that("subset_h5ad works with dense matrix", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "array")

  cell_idx <- c(10L, 5L, 50L, 75L, 1L)
  gene_idx <- c(3L, 1L, 20L, 40L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = gene_idx,
    .chunk.size  = 2L,
    .verbose     = FALSE
  )

  expect_true(file.exists(dst))

  set.seed(42)
  original_dense <- matrix(rpois(100 * 50, lambda = 0.3), nrow = 100, ncol = 50)
  expected_sub <- original_dense[cell_idx, gene_idx]

  result_dense <- .read_output_matrix_as_dense(dst, length(cell_idx),
                                               length(gene_idx))
  expect_equal(result_dense, expected_sub)
})

test_that("subset_h5ad works with CSC matrix", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csc_matrix")

  cell_idx <- c(10L, 5L, 50L, 75L, 1L)
  gene_idx <- c(3L, 1L, 20L, 40L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = gene_idx,
    .chunk.size  = 2L,
    .verbose     = FALSE
  )

  expect_true(file.exists(dst))

  set.seed(42)
  original_dense <- matrix(rpois(100 * 50, lambda = 0.3), nrow = 100, ncol = 50)
  expected_sub <- original_dense[cell_idx, gene_idx]

  result_dense <- .read_output_matrix_as_dense(dst, length(cell_idx),
                                               length(gene_idx))
  expect_equal(result_dense, expected_sub)
})

test_that("subset_h5ad preserves obs metadata correctly", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  cell_idx <- c(10L, 5L, 50L, 75L, 1L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .verbose     = FALSE
  )

  # Cell names
  obs_index <- rhdf5::h5read(dst, "/obs/_index")
  src_cell_names <- paste0("cell_", seq_len(100))
  expect_equal(c(obs_index), src_cell_names[cell_idx])

  # sample_id column
  out_sample <- rhdf5::h5read(dst, "/obs/sample_id")
  src_sample <- paste0("sample_", rep(1:5, length.out = 100))
  expect_equal(c(out_sample), src_sample[cell_idx])

  # Categorical cell_type column
  out_codes <- rhdf5::h5read(dst, "/obs/cell_type/codes")
  set.seed(123)
  src_codes <- sample(0:2, 100, replace = TRUE)
  expect_equal(c(out_codes), src_codes[cell_idx])

  # Categories should be copied in full
  out_cats <- rhdf5::h5read(dst, "/obs/cell_type/categories")
  expect_equal(c(out_cats), c("T_cell", "B_cell", "Monocyte"))

  # Encoding attributes
  obs_attrs <- rhdf5::h5readAttributes(dst, "/obs")
  expect_equal(obs_attrs[["encoding-type"]], "dataframe")
  expect_equal(obs_attrs[["_index"]], "_index")
})

test_that("subset_h5ad preserves var metadata correctly", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  gene_idx <- c(40L, 10L, 1L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = seq_len(100L),
    .gene.idx    = gene_idx,
    .verbose     = FALSE
  )

  var_index <- rhdf5::h5read(dst, "/var/_index")
  src_gene_names <- paste0("gene_", seq_len(50))
  expect_equal(c(var_index), src_gene_names[gene_idx])

  var_attrs <- rhdf5::h5readAttributes(dst, "/var")
  expect_equal(var_attrs[["encoding-type"]], "dataframe")
  expect_equal(var_attrs[["_index"]], "_index")
})

test_that("subset_h5ad copies obsm embeddings", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix",
                   with_obsm = TRUE)

  cell_idx <- c(10L, 5L, 50L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .copy.obsm   = TRUE,
    .verbose     = FALSE
  )

  # Read source PCA
  set.seed(99)
  src_pca <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)

  out_pca <- rhdf5::h5read(dst, "/obsm/X_pca")
  expect_equal(nrow(out_pca), length(cell_idx))
  expect_equal(ncol(out_pca), 10L)
  expect_equal(out_pca, src_pca[cell_idx, , drop = FALSE])
})

test_that("subset_h5ad gene_idx = NULL keeps all genes", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  cell_idx <- c(1L, 10L, 20L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = NULL,
    .verbose     = FALSE
  )

  var_index <- rhdf5::h5read(dst, "/var/_index")
  expect_equal(length(var_index), 50L)

  set.seed(42)
  original_dense <- matrix(rpois(100 * 50, lambda = 0.3), nrow = 100, ncol = 50)
  expected_sub <- original_dense[cell_idx, ]

  result_dense <- .read_output_matrix_as_dense(dst, length(cell_idx), 50L)
  expect_equal(result_dense, expected_sub)
})

test_that("subset_h5ad handles layers", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix",
                   n_layers = 1L)

  cell_idx <- c(1L, 10L, 20L)
  gene_idx <- c(5L, 15L, 25L)

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = gene_idx,
    .verbose     = FALSE
  )

  # Check layer was copied
  fid <- rhdf5::H5Fopen(dst)
  has_layer <- rhdf5::H5Lexists(fid, "/layers/layer_1")
  rhdf5::H5Fclose(fid)
  expect_true(has_layer)

  # Verify layer values
  set.seed(201)
  src_layer <- matrix(runif(100 * 50), nrow = 100, ncol = 50)
  expected_layer <- src_layer[cell_idx, gene_idx]

  out_layer <- rhdf5::h5read(dst, "/layers/layer_1")
  expect_equal(out_layer, expected_layer, tolerance = 1e-10)
})

test_that("subset_h5ad validates inputs", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  # Non-existent source

  expect_error(
    subset_h5ad(.source.path = "nonexistent.h5ad",
                .dest.path   = dst,
                .cell.idx    = 1:5),
    regexp = NULL
  )

  # Same source and dest
  expect_error(
    subset_h5ad(.source.path = src,
                .dest.path   = src,
                .cell.idx    = 1:5),
    regexp = NULL
  )

  # Empty cell_idx
  expect_error(
    subset_h5ad(.source.path = src,
                .dest.path   = dst,
                .cell.idx    = integer(0)),
    regexp = NULL
  )

  # Out-of-range cell_idx
  expect_error(
    subset_h5ad(.source.path = src,
                .dest.path   = dst,
                .cell.idx    = c(1L, 200L)),
    regexp = NULL
  )
})

test_that("subset_h5ad copies /uns metadata", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = 1:10,
    .copy.uns    = TRUE,
    .verbose     = FALSE
  )

  out_dn <- rhdf5::h5read(dst, "/uns/dataset_name")
  expect_equal(c(out_dn), "test_dataset")
})

test_that("subset_h5ad output is readable by anndataR", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("anndataR")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = c(1L, 10L, 20L, 30L),
    .gene.idx    = c(1L, 5L, 10L),
    .verbose     = FALSE
  )

  adata <- anndataR::read_h5ad(dst, as = "HDF5AnnData")
  expect_equal(nrow(adata), 4L)
  expect_equal(ncol(adata), 3L)
})

test_that("subset_h5ad CSR with gene_idx=NULL and reordered cells", {
  skip_if_not_installed("rhdf5")

  src <- tempfile(fileext = ".h5ad")
  dst <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(src, dst)), add = TRUE)

  create_test_h5ad(src, n_obs = 100L, n_var = 50L, encoding = "csr_matrix")

  # Reverse order of first 20 cells
  cell_idx <- 20:1

  subset_h5ad(
    .source.path = src,
    .dest.path   = dst,
    .cell.idx    = cell_idx,
    .gene.idx    = NULL,
    .chunk.size  = 3L,
    .verbose     = FALSE
  )

  set.seed(42)
  original_dense <- matrix(rpois(100 * 50, lambda = 0.3), nrow = 100, ncol = 50)
  expected_sub <- original_dense[cell_idx, ]

  result_dense <- .read_output_matrix_as_dense(dst, length(cell_idx), 50L)
  expect_equal(result_dense, expected_sub)
})
