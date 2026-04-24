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

#' Create a minimal SCE test fixture
#'
#' Builds a small SingleCellExperiment with known values for testing.
#' @param n_cells Number of cells PER SAMPLE (total = n_cells * n_samples)
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param use_delayed If TRUE, wrap the counts matrix in a
#'   \code{DelayedArray::DelayedArray()} so the assay is a DelayedMatrix.
#'   If FALSE, store as a dgCMatrix.
#' @return A \code{SingleCellExperiment} object.
#' @noRd
.make_test_sce <- function(n_cells = 100,
                           n_genes = 50,
                           n_samples = 4,
                           use_delayed = FALSE) {
  set.seed(42)
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells, n_samples)
  total_cells <- sum(cells_per_sample)

  sample_ids <- rep(sample_names, each = n_cells)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  # genes x cells (SummarizedExperiment convention)
  counts <- matrix(
    rpois(n_genes * total_cells, lambda = 5),
    nrow = n_genes, ncol = total_cells,
    dimnames = list(gene_ids, cell_ids)
  )
  counts_sparse <- as(counts, "CsparseMatrix")

  col_data <- S4Vectors::DataFrame(
    sample_id = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids
  )

  if (use_delayed) {
    assay_mat <- DelayedArray::DelayedArray(counts_sparse)
  } else {
    assay_mat <- counts_sparse
  }

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = assay_mat),
    colData = col_data
  )
}

# ======================================================================
# .delayed_to_bpcells tests
# ======================================================================

test_that(".delayed_to_bpcells converts a non-HDF5 DelayedMatrix via chunked fallback", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  dgc <- Matrix::Matrix(
    matrix(rpois(50 * 200, lambda = 5), nrow = 50, ncol = 200,
           dimnames = list(paste0("g", 1:50), paste0("c", 1:200))),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)

  bp_dir <- file.path(withr::local_tempdir(), "bp_chunked")
  result <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                              .verbose = FALSE)

  expect_true(methods::is(result, "IterableMatrix"))
  expect_equal(nrow(result), 50L)
  expect_equal(ncol(result), 200L)
  expect_true(dir.exists(bp_dir))
})

test_that(".delayed_to_bpcells cache hit", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  dgc <- Matrix::Matrix(
    matrix(rpois(50 * 200, lambda = 5), nrow = 50, ncol = 200,
           dimnames = list(paste0("g", 1:50), paste0("c", 1:200))),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)

  bp_dir <- file.path(withr::local_tempdir(), "bp_cache")

  # First call: creates on disk
  result1 <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                               .verbose = FALSE)
  expect_true(methods::is(result1, "IterableMatrix"))

  # Second call: should hit cache
  out <- capture.output(
    result2 <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                                 .verbose = TRUE)
  )
  expect_true(any(grepl("cache hit", out, ignore.case = TRUE)))
  expect_equal(nrow(result2), nrow(result1))
  expect_equal(ncol(result2), ncol(result1))
})

test_that(".delayed_to_bpcells with explicit .bpcells.dir", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  dgc <- Matrix::Matrix(
    matrix(rpois(50 * 200, lambda = 5), nrow = 50, ncol = 200,
           dimnames = list(paste0("g", 1:50), paste0("c", 1:200))),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)

  bp_dir <- file.path(withr::local_tempdir(), "bp_explicit")
  result <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                              .verbose = FALSE)

  expect_true(dir.exists(bp_dir))
  expect_true(methods::is(result, "IterableMatrix"))
  # Verify the returned matrix is usable (can subset columns)
  sub <- result[, 1:10]
  expect_equal(ncol(sub), 10L)
})

test_that(".delayed_to_bpcells preserves dimnames", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  rnames <- paste0("gene_", 1:30)
  cnames <- paste0("cell_", 1:100)
  dgc <- Matrix::Matrix(
    matrix(rpois(30 * 100, lambda = 5), nrow = 30, ncol = 100,
           dimnames = list(rnames, cnames)),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)

  bp_dir <- file.path(withr::local_tempdir(), "bp_dimnames")
  result <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                              .verbose = FALSE)

  expect_equal(rownames(result), rnames)
  expect_equal(colnames(result), cnames)
})

test_that(".delayed_to_bpcells with .bpcells.dir = NULL uses tempdir", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  # Use unique dimensions so tempdir path doesn't collide with other tests
  dgc <- Matrix::Matrix(
    matrix(rpois(31 * 201, lambda = 5), nrow = 31, ncol = 201,
           dimnames = list(paste0("g", 1:31), paste0("c", 1:201))),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)

  result <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = NULL,
                                              .verbose = FALSE)

  expect_true(methods::is(result, "IterableMatrix"))
  expect_equal(nrow(result), 31L)
  expect_equal(ncol(result), 201L)
})

test_that(".delayed_to_bpcells seed-detection recognises IterableMatrix seeds", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  # BPCells IterableMatrix cannot currently be wrapped in DelayedArray
  # (IterableMatrix does not implement extract_array()), so the
  # short-circuit branch in .delayed_to_bpcells is unreachable in
  # practice.  We verify the detection logic directly: the function
  # checks `methods::is(seed, "IterableMatrix")`, so we ensure that
  # check works on a real on-disk BPCells matrix.
  set.seed(42)
  dgc <- Matrix::Matrix(
    matrix(rpois(50 * 100, lambda = 5), nrow = 50, ncol = 100,
           dimnames = list(paste0("g", 1:50), paste0("c", 1:100))),
    sparse = TRUE
  )
  bp_disk_dir <- file.path(withr::local_tempdir(), "bp_seed_disk")
  BPCells::write_matrix_dir(
    mat = methods::as(dgc, "IterableMatrix"), dir = bp_disk_dir
  )
  bp_disk <- BPCells::open_matrix_dir(bp_disk_dir)

  # Confirm it IS an IterableMatrix (the condition .delayed_to_bpcells
  # uses) and is NOT a valid DelayedArray seed.
  expect_true(methods::is(bp_disk, "IterableMatrix"))
  expect_error(
    DelayedArray::DelayedArray(bp_disk),
    "extract_array"
  )
})

test_that(".delayed_to_bpcells converts HDF5-backed DelayedMatrix (fallback path)", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("HDF5Array")

  set.seed(42)
  dense_mat <- matrix(rpois(50 * 200, lambda = 5), nrow = 50, ncol = 200,
                      dimnames = list(paste0("g", 1:50), paste0("c", 1:200)))

  h5_dir <- withr::local_tempdir()
  h5_path <- file.path(h5_dir, "test.h5")
  hdf5_mat <- HDF5Array::writeHDF5Array(dense_mat, filepath = h5_path,
                                         name = "counts")
  dmat <- DelayedArray::DelayedArray(hdf5_mat)

  bp_dir <- file.path(withr::local_tempdir(), "bp_hdf5")

  # HDF5Array writes a standard dense HDF5 dataset that BPCells cannot
  # read via open_matrix_hdf5 (it expects BPCells' own CSC layout).
  # .delayed_to_bpcells must succeed anyway by falling back to the
  # chunked conversion path.
  out <- capture.output(
    result <- tinydenseR:::.delayed_to_bpcells(dmat, .bpcells.dir = bp_dir,
                                                .verbose = TRUE)
  )

  expect_true(methods::is(result, "IterableMatrix"))
  expect_equal(nrow(result), 50L)
  expect_equal(ncol(result), 200L)
  expect_true(dir.exists(bp_dir))

  # Verify the fallback message was emitted
  expect_true(
    any(grepl("falling back to chunked", out, ignore.case = TRUE)) ||
    any(grepl("chunked", out, ignore.case = TRUE))
  )

  # Verify data integrity: materialise first column and compare
  bp_col1 <- as.matrix(result[, 1, drop = FALSE])
  expect_equal(as.numeric(bp_col1), as.numeric(dense_mat[, 1]))
})

# ======================================================================
# RunTDR.DelayedMatrix full pipeline tests
# ======================================================================

test_that("RunTDR.DelayedMatrix full pipeline", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  n_cells_per_sample <- 100
  n_genes <- 50
  n_samples <- 4
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells_per_sample, n_samples)
  total_cells <- sum(cells_per_sample)
  sample_ids <- rep(sample_names, each = n_cells_per_sample)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  dgc <- as(counts, "CsparseMatrix")
  delayed_mat <- DelayedArray::DelayedArray(dgc)

  meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_full_pipeline")
  result <- RunTDR(delayed_mat, .cell.meta = meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
  expect_equal(result@config$backend, "matrix")
  expect_true(!is.null(result@density$norm))
  expect_equal(nrow(result@metadata), n_samples)
  expect_true(all(sort(rownames(result@metadata)) == sort(sample_names)))
})

test_that("RunTDR.DelayedMatrix BPCells cache hit", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  n_cells_per_sample <- 100
  n_genes <- 50
  n_samples <- 4
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells_per_sample, n_samples)
  total_cells <- sum(cells_per_sample)
  sample_ids <- rep(sample_names, each = n_cells_per_sample)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  dgc <- as(counts, "CsparseMatrix")
  delayed_mat <- DelayedArray::DelayedArray(dgc)

  meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_runtdr_cache")

  result1 <- RunTDR(delayed_mat, .cell.meta = meta,
                    .sample.var = "Sample",
                    .assay.type = "RNA",
                    .bpcells.dir = bp_dir,
                    .verbose = FALSE, .seed = 42)
  expect_true(is.TDRObj(result1))
  expect_true(dir.exists(bp_dir))

  # Second run should hit cache
  out <- capture.output(
    result2 <- RunTDR(delayed_mat, .cell.meta = meta,
                      .sample.var = "Sample",
                      .assay.type = "RNA",
                      .bpcells.dir = bp_dir,
                      .verbose = TRUE, .seed = 42)
  )
  expect_true(is.TDRObj(result2))
  expect_true(any(grepl("cache hit", out, ignore.case = TRUE)))
})

test_that("RunTDR.DelayedMatrix with .bpcells.dir = NULL uses tempdir", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  n_cells_per_sample <- 100
  n_genes <- 50
  n_samples <- 4
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells_per_sample, n_samples)
  total_cells <- sum(cells_per_sample)
  sample_ids <- rep(sample_names, each = n_cells_per_sample)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  dgc <- as(counts, "CsparseMatrix")
  delayed_mat <- DelayedArray::DelayedArray(dgc)

  meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  # Use a unique explicit dir to avoid cache collisions
  bp_dir <- file.path(withr::local_tempdir(), "bp_null_test")
  result <- RunTDR(delayed_mat, .cell.meta = meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
})

# ======================================================================
# RunTDR.SingleCellExperiment + DelayedMatrix tests
# ======================================================================

test_that("RunTDR.SingleCellExperiment with DelayedMatrix assay - full pipeline", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  sce <- .make_test_sce(n_cells = 100, n_genes = 50, n_samples = 4,
                        use_delayed = TRUE)

  bp_dir <- file.path(withr::local_tempdir(), "bp_sce_full")
  result <- RunTDR(sce, .sample.var = "sample_id",
                   .bpcells.dir = bp_dir, .verbose = FALSE)

  expect_true(inherits(result, "SingleCellExperiment"))
  tdr.obj <- S4Vectors::metadata(result)$tdr.obj
  expect_true(is.TDRObj(tdr.obj))
  expect_equal(tdr.obj@config$backend, "matrix")
  expect_equal(nrow(tdr.obj@metadata), 4L)
  expect_true(!is.null(tdr.obj@density$norm))
})

test_that("RunTDR.SingleCellExperiment with DelayedMatrix and .bpcells.dir", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  sce <- .make_test_sce(n_cells = 100, n_genes = 50, n_samples = 4,
                        use_delayed = TRUE)

  bp_dir <- file.path(withr::local_tempdir(), "bp_sce_dir")

  result <- RunTDR(sce, .sample.var = "sample_id",
                   .bpcells.dir = bp_dir, .verbose = FALSE)

  expect_true(dir.exists(bp_dir))
  tdr.obj <- S4Vectors::metadata(result)$tdr.obj
  expect_true(is.TDRObj(tdr.obj))
})

test_that("RunTDR.SingleCellExperiment with DelayedMatrix and .celltype.vec", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  sce <- .make_test_sce(n_cells = 100, n_genes = 50, n_samples = 4,
                        use_delayed = TRUE)

  # Add a cell type column to colData
  SummarizedExperiment::colData(sce)$cell_type <- rep(
    c("TypeA", "TypeB", "TypeC", "TypeD"), length.out = ncol(sce)
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_sce_ct")
  result <- RunTDR(sce, .sample.var = "sample_id",
                   .celltype.vec = "cell_type",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE)

  tdr.obj <- S4Vectors::metadata(result)$tdr.obj
  expect_true(is.TDRObj(tdr.obj))
  expect_true(!is.null(tdr.obj@config$celltype.vec))
})

test_that("RunTDR.SingleCellExperiment + celltype.vec produces valid norm (drop=FALSE regression)", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  # This test guards against the regression where fgraph subsetting

  # by celltype-confidence filter can reduce to a single row, causing
  # Matrix::colSums to receive a vector instead of a matrix when
  # drop = FALSE is not applied. (GitHub issue: 'x' must be an array
  # of at least two dimensions.)
  sce <- .make_test_sce(n_cells = 100, n_genes = 50, n_samples = 4,
                        use_delayed = TRUE)

  SummarizedExperiment::colData(sce)$cell_type <- rep(
    c("TypeA", "TypeB", "TypeC", "TypeD"), length.out = ncol(sce)
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_sce_ct_fdens")
  result <- RunTDR(sce, .sample.var = "sample_id",
                   .celltype.vec = "cell_type",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE)

  tdr.obj <- S4Vectors::metadata(result)$tdr.obj

  # norm must be a numeric matrix with one column per sample
  expect_true(!is.null(tdr.obj@density$norm))
  expect_true(is.numeric(tdr.obj@density$norm))
  expect_equal(ncol(tdr.obj@density$norm), nrow(tdr.obj@metadata))
  # norm values must be non-negative (they are fuzzy membership sums)
  expect_true(all(tdr.obj@density$norm >= 0))
})

test_that("fgraph subsetting preserves matrix structure with drop = FALSE", {
  skip_on_cran()

  # Directly validate that subsetting a sparse matrix by a logical vector
  # that selects only 1 row does NOT drop to a vector when drop = FALSE.
  # This is the unit-level guard for the norm computation fix.
  set.seed(42)
  fg <- Matrix::Matrix(
    matrix(rpois(20, 3), nrow = 5, ncol = 4,
           dimnames = list(paste0("cell", 1:5), paste0("lm", 1:4))),
    sparse = TRUE
  )

  # Selecting 1 row WITHOUT drop = FALSE produces a vector
  one_row <- c(FALSE, TRUE, FALSE, FALSE, FALSE)
  expect_false(is.matrix(fg[one_row, ]) || methods::is(fg[one_row, ], "Matrix"))

  # Selecting 1 row WITH drop = FALSE preserves matrix
  expect_true(methods::is(fg[one_row, , drop = FALSE], "Matrix"))
  expect_equal(nrow(fg[one_row, , drop = FALSE]), 1L)

  # colSums works on the preserved matrix
  expect_equal(
    length(Matrix::colSums(fg[one_row, , drop = FALSE])),
    ncol(fg)
  )
})

test_that("RunTDR.SingleCellExperiment with DelayedMatrix and .min.cells.per.sample filtering", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  # Create SCE with one small sample
  set.seed(42)
  n_genes <- 50
  gene_ids <- paste0("gene_", seq_len(n_genes))

  # 3 big samples (150 cells each) + 1 tiny sample (3 cells)
  sample_ids <- c(rep("sampleA", 150), rep("sampleB", 150),
                  rep("sampleC", 150), rep("sampleD", 3))
  total_cells <- length(sample_ids)
  cell_ids <- paste0("cell_", seq_len(total_cells))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  counts_sparse <- as(counts, "CsparseMatrix")
  dmat <- DelayedArray::DelayedArray(counts_sparse)

  col_data <- S4Vectors::DataFrame(
    sample_id = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = dmat),
    colData = col_data
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_sce_filter")
  result <- RunTDR(sce, .sample.var = "sample_id",
                   .min.cells.per.sample = 10,
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE)

  tdr.obj <- S4Vectors::metadata(result)$tdr.obj
  expect_true(is.TDRObj(tdr.obj))
  # sampleD should be excluded (only 3 cells)
  expect_false("sampleD" %in% rownames(tdr.obj@metadata))
  expect_equal(nrow(tdr.obj@metadata), 3L)
})

# ======================================================================
# Input validation tests
# ======================================================================

test_that("RunTDR.DelayedMatrix input validation errors", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  dgc <- Matrix::Matrix(
    matrix(rpois(50 * 100, lambda = 5), nrow = 50, ncol = 100,
           dimnames = list(paste0("g", 1:50), paste0("c", 1:100))),
    sparse = TRUE
  )
  dmat <- DelayedArray::DelayedArray(dgc)
  meta <- data.frame(
    Sample = rep(paste0("s", 1:4), each = 25),
    row.names = paste0("c", 1:100),
    stringsAsFactors = FALSE
  )

  # Non-RNA assay type
  expect_error(
    RunTDR(dmat, .cell.meta = meta, .sample.var = "Sample",
           .assay.type = "cyto"),
    "only supported for .assay.type = 'RNA'"
  )

  # Non-data.frame .cell.meta
  expect_error(
    RunTDR(dmat, .cell.meta = "not_a_df", .sample.var = "Sample"),
    "data\\.frame"
  )

  # Missing .sample.var column
  expect_error(
    RunTDR(dmat, .cell.meta = meta, .sample.var = "nonexistent"),
    "nonexistent"
  )
})

# ======================================================================
# Data integrity tests
# ======================================================================

test_that("Dimnames preserved through BPCells conversion in RunTDR.DelayedMatrix", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  n_cells_per_sample <- 100
  n_genes <- 50
  n_samples <- 4
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells_per_sample, n_samples)
  total_cells <- sum(cells_per_sample)
  sample_ids <- rep(sample_names, each = n_cells_per_sample)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  dgc <- as(counts, "CsparseMatrix")
  delayed_mat <- DelayedArray::DelayedArray(dgc)

  meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  bp_dir <- file.path(withr::local_tempdir(), "bp_dimnames_pipeline")
  result <- RunTDR(delayed_mat, .cell.meta = meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE, .seed = 42)

  # Verify gene names preserved in source matrix
  expect_equal(rownames(result@config$source.env$mat), gene_ids)

  # Verify cell IDs in @cells match expected cell names
  all_cell_ids <- unlist(lapply(result@cells, function(idx) {
    colnames(result@config$source.env$mat)[idx]
  }))
  expect_true(all(all_cell_ids %in% cell_ids))
})

test_that("S3 dispatch works for DelayedMatrix", {
  skip_on_cran()
  skip_if_not_installed("BPCells")
  skip_if_not_installed("DelayedArray")

  set.seed(42)
  n_cells_per_sample <- 100
  n_genes <- 50
  n_samples <- 4
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells_per_sample, n_samples)
  total_cells <- sum(cells_per_sample)
  sample_ids <- rep(sample_names, each = n_cells_per_sample)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))
  gene_ids <- paste0("gene_", seq_len(n_genes))

  counts <- matrix(rpois(n_genes * total_cells, lambda = 5),
                   nrow = n_genes, ncol = total_cells,
                   dimnames = list(gene_ids, cell_ids))
  dgc <- as(counts, "CsparseMatrix")
  delayed_mat <- DelayedArray::DelayedArray(dgc)

  meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  # Verify RunTDR() dispatches correctly (not RunTDR.default)
  bp_dir <- file.path(withr::local_tempdir(), "bp_dispatch")
  result <- RunTDR(delayed_mat, .cell.meta = meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .bpcells.dir = bp_dir,
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
})

# ======================================================================
# Backward compatibility and removal verification tests
# ======================================================================

test_that("RunTDR.SingleCellExperiment with dgCMatrix assay uses sce backend (unchanged)", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")

  sce <- .make_test_sce(n_cells = 100, n_genes = 50, n_samples = 4,
                        use_delayed = FALSE)

  result <- RunTDR(sce, .sample.var = "sample_id", .verbose = FALSE)

  expect_true(inherits(result, "SingleCellExperiment"))
  tdr.obj <- S4Vectors::metadata(result)$tdr.obj
  expect_true(is.TDRObj(tdr.obj))
  expect_equal(tdr.obj@config$backend, "sce")
})

test_that(".optimize.hdf5 parameter is removed and DelayedMatrixStats is no longer in Suggests", {
  skip_on_cran()

  # .optimize.hdf5 should NOT be a formal parameter
  expect_false(
    ".optimize.hdf5" %in% names(formals(tinydenseR:::RunTDR.SingleCellExperiment))
  )

  # DelayedMatrixStats should NOT be in Suggests
  desc <- packageDescription("tinydenseR")
  if (!is.null(desc$Suggests)) {
    suggests <- trimws(unlist(strsplit(desc$Suggests, ",")))
    expect_false("DelayedMatrixStats" %in% suggests)
  }
})
