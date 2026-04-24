#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ======================================================================
# Helper: create a dgCMatrix + cell.meta for testing
# ======================================================================

.make_matrix_test_data <- function(n_cells = 100, n_markers = 50,
                                   n_samples = 4, assay.type = "RNA") {
  sample_names <- paste0("sample", seq_len(n_samples))
  cells_per_sample <- rep(n_cells, n_samples)
  total_cells <- sum(cells_per_sample)
  marker_names <- paste0("marker_", seq_len(n_markers))

  sample_ids <- rep(sample_names, each = n_cells)
  cell_ids <- paste0(sample_ids, "_cell_", sequence(cells_per_sample))

  if (assay.type == "cyto") {
    # cells x markers
    mat <- Matrix::Matrix(
      data = matrix(runif(total_cells * n_markers),
                    nrow = total_cells, ncol = n_markers,
                    dimnames = list(cell_ids, marker_names)),
      sparse = TRUE
    )
  } else {
    # features x cells
    mat <- Matrix::Matrix(
      data = matrix(rpois(n_markers * total_cells, lambda = 5),
                    nrow = n_markers, ncol = total_cells,
                    dimnames = list(marker_names, cell_ids)),
      sparse = TRUE
    )
  }

  cell.meta <- data.frame(
    Sample = sample_ids,
    group = rep(c("A", "B"), length.out = total_cells),
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )

  list(mat = mat, cell.meta = cell.meta)
}

# ======================================================================
# S3 dispatch tests
# ======================================================================

test_that("RunTDR dispatches to dgCMatrix method", {
  expect_true(is.function(getS3method("RunTDR", "dgCMatrix")))
})

test_that("RunTDR dispatches to DelayedMatrix method", {
  expect_true(is.function(getS3method("RunTDR", "DelayedMatrix")))
})

test_that("RunTDR dispatches to IterableMatrix method", {
  expect_true(is.function(getS3method("RunTDR", "IterableMatrix")))
})

# ======================================================================
# RunTDR.dgCMatrix tests
# ======================================================================

test_that("RunTDR.dgCMatrix rejects cyto assay type", {
  td <- .make_matrix_test_data(n_cells = 50, n_markers = 5,
                               n_samples = 2, assay.type = "cyto")
  expect_error(RunTDR(td$mat, .cell.meta = td$cell.meta,
                      .sample.var = "Sample",
                      .markers = paste0("marker_", 1:5),
                      .assay.type = "cyto"),
               "only supported for .assay.type = 'RNA'")
})

test_that("RunTDR.dgCMatrix runs full pipeline (RNA)", {
  skip_on_cran()
  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  result <- RunTDR(td$mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
  expect_equal(result@config$backend, "matrix")
  expect_true(!is.null(result@density$norm))
  expect_equal(nrow(result@metadata), 4L)
})

test_that("RunTDR.dgCMatrix uses index-based .cells (not file paths)", {
  skip_on_cran()
  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  result <- RunTDR(td$mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  # .cells should be integer index vectors, not file paths
  expect_true(is.integer(result@cells[[1]]))
  expect_true(all(vapply(result@cells, is.integer, logical(1))))
})

test_that("RunTDR.dgCMatrix stores a locked source env", {
  skip_on_cran()
  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  result <- RunTDR(td$mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  env <- result@config$source.env
  expect_true(is.environment(env))
  expect_true(environmentIsLocked(env) || bindingIsLocked("mat", env))
  expect_error(env$mat <- "overwrite", "cannot change")
})

test_that(".get_sample_matrix works with matrix backend (RNA)", {
  skip_on_cran()
  td <- .make_matrix_test_data(n_cells = 350, n_markers = 50,
                               n_samples = 3, assay.type = "RNA")

  result <- RunTDR(td$mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  sample1_mat <- tinydenseR:::.get_sample_matrix(NULL, result, 1)
  expect_equal(ncol(sample1_mat), lengths(result@cells)[[1]])
  expect_equal(nrow(sample1_mat), 50L)
})

# ======================================================================
# RunTDR.DelayedMatrix tests
# ======================================================================

test_that("RunTDR.DelayedMatrix runs full pipeline", {
  skip_on_cran()
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("BPCells")

  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  delayed_mat <- DelayedArray::DelayedArray(as.matrix(td$mat))

  result <- RunTDR(delayed_mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
  expect_equal(result@config$backend, "matrix")
  expect_true(!is.null(result@density$norm))
})

# ======================================================================
# RunTDR.IterableMatrix tests
# ======================================================================

test_that("RunTDR.IterableMatrix runs full pipeline", {
  skip_on_cran()
  skip_if_not_installed("BPCells")

  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  bp_mat <- methods::as(td$mat, "IterableMatrix")

  result <- RunTDR(bp_mat, .cell.meta = td$cell.meta,
                   .sample.var = "Sample",
                   .assay.type = "RNA",
                   .verbose = FALSE, .seed = 42)

  expect_true(is.TDRObj(result))
  expect_equal(result@config$backend, "matrix")
  expect_true(!is.null(result@density$norm))
})

# ======================================================================
# Input validation error tests
# ======================================================================

test_that("RunTDR.dgCMatrix errors on non-data.frame .cell.meta", {
  mat <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1:3,
                              dimnames = list(paste0("g", 1:3),
                                              paste0("c", 1:3)))
  expect_error(RunTDR(mat, .cell.meta = "not_a_df",
                      .sample.var = "x",
                      .assay.type = "RNA"),
               "\\.cell\\.meta must be a data\\.frame")
})

test_that("RunTDR.dgCMatrix errors on missing .sample.var column", {
  td <- .make_matrix_test_data(n_cells = 10, n_markers = 3,
                               n_samples = 2, assay.type = "RNA")
  expect_error(RunTDR(td$mat, .cell.meta = td$cell.meta,
                      .sample.var = "NonExistent",
                      .assay.type = "RNA"),
               "not found in \\.cell\\.meta")
})

test_that("RunTDR.dgCMatrix errors on missing matrix dimnames", {
  mat <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1:3,
                              dims = c(3, 3))
  cell.meta <- data.frame(Sample = c("s1", "s1", "s2"),
                           row.names = paste0("c", 1:3))
  expect_error(RunTDR(mat, .cell.meta = cell.meta,
                      .sample.var = "Sample",
                      .assay.type = "RNA"),
               "must have colnames")
})

test_that("RunTDR.dgCMatrix errors when no cell overlap", {
  mat <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1:3,
                              dimnames = list(paste0("g", 1:3),
                                              paste0("c", 1:3)))
  cell.meta <- data.frame(Sample = c("s1", "s2"),
                           row.names = c("x1", "x2"))
  expect_error(RunTDR(mat, .cell.meta = cell.meta,
                      .sample.var = "Sample",
                      .assay.type = "RNA"),
               "No overlap")
})

test_that("RunTDR.dgCMatrix errors when all samples below min cells", {
  td <- .make_matrix_test_data(n_cells = 5, n_markers = 3,
                               n_samples = 2, assay.type = "RNA")
  expect_error(RunTDR(td$mat, .cell.meta = td$cell.meta,
                      .sample.var = "Sample",
                      .assay.type = "RNA",
                      .min.cells.per.sample = 1000),
               "No samples have")
})

# ======================================================================
# Equivalence: matrix backend vs files backend (via setup.tdr.obj)
# ======================================================================

test_that("matrix backend produces same cell counts as files backend", {
  skip_on_cran()
  td <- .make_matrix_test_data(n_cells = 100, n_markers = 50,
                               n_samples = 4, assay.type = "RNA")

  result_mat <- RunTDR(td$mat, .cell.meta = td$cell.meta,
                       .sample.var = "Sample",
                       .assay.type = "RNA",
                       .verbose = FALSE, .seed = 42)

  # Files-based equivalent
  sample_names <- paste0("sample", 1:4)
  mat_list <- lapply(stats::setNames(sample_names, sample_names), function(s) {
    idx <- which(td$cell.meta$Sample == s)
    td$mat[, idx, drop = FALSE]
  })
  files <- lapply(mat_list, function(m) {
    f <- tempfile(fileext = ".RDS")
    saveRDS(m, f, compress = FALSE)
    f
  })
  on.exit(lapply(files, unlink), add = TRUE)

  sample_meta <- data.frame(
    Sample = sample_names,
    group = c("A", "B", "A", "B"),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
  tdr_files <- setup.tdr.obj(.cells = files, .meta = sample_meta,
                             .assay.type = "RNA", .verbose = FALSE)
  tdr_files <- RunTDR(tdr_files, .verbose = FALSE, .seed = 42)

  expect_equal(result_mat@config$sampling$n.cells,
               tdr_files@config$sampling$n.cells)
})
