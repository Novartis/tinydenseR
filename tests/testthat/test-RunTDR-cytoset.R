#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ======================================================================
# Helper: create a cytoset + pdata for testing
# ======================================================================

.make_cytoset_test_data <- function(n_cells = 100, n_markers = 5,
                                    n_samples = 4) {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  sample_names <- paste0("sample", seq_len(n_samples))
  tmp_dir <- tempfile("cytoset_test_")
  dir.create(tmp_dir)

  fcs_paths <- vapply(sample_names, function(sn) {
    mat <- matrix(runif(n_cells * n_markers),
                  nrow = n_cells, ncol = n_markers,
                  dimnames = list(paste0("event_", seq_len(n_cells)),
                                  paste0("Marker", seq_len(n_markers))))
    ff <- flowCore::flowFrame(exprs = mat)
    fpath <- file.path(tmp_dir, paste0(sn, ".fcs"))
    flowCore::write.FCS(ff, filename = fpath)
    fpath
  }, character(1))

  cs <- flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(fcs_paths, sample_names))

  pd <- flowWorkspace::pData(cs)
  pd$Sample <- sample_names
  pd$group <- rep(c("A", "B"), length.out = n_samples)
  flowWorkspace::pData(cs) <- pd

  list(cs = cs, pdata = pd, tmp_dir = tmp_dir)
}

# ======================================================================
# S3 dispatch
# ======================================================================

test_that("RunTDR.cytoset S3 method exists", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")
  expect_true(is.function(getS3method("RunTDR", "cytoset")))
})

# ======================================================================
# Input validation
# ======================================================================

test_that("RunTDR.cytoset errors when .sample.var not in pData", {
  td <- .make_cytoset_test_data()
  expect_error(
    RunTDR(td$cs, .sample.var = "nonexistent"),
    "not found in pData"
  )
})

test_that("RunTDR.cytoset errors when .assay.type is not cyto", {
  td <- .make_cytoset_test_data()
  expect_error(
    RunTDR(td$cs, .sample.var = "Sample", .assay.type = "RNA"),
    "must be 'cyto'"
  )
})

test_that("RunTDR.cytoset errors when .markers not in channels", {
  td <- .make_cytoset_test_data()
  expect_error(
    RunTDR(td$cs, .sample.var = "Sample", .markers = c("FakeMarker")),
    "not found in channels"
  )
})

# ======================================================================
# Full pipeline
# ======================================================================

test_that("RunTDR.cytoset full pipeline returns TDRObj", {
  skip_on_cran()
  td <- .make_cytoset_test_data()

  tdr <- RunTDR(td$cs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .verbose = FALSE,
                .seed = 42)

  expect_true(is.TDRObj(tdr))
  expect_equal(tdr@config$backend, "cyto")
  expect_true(length(tdr@density$norm) > 0)
  expect_equal(nrow(tdr@metadata), length(tdr@cells))
})

# ======================================================================
# Locked environment
# ======================================================================

test_that("RunTDR.cytoset locks source env", {
  skip_on_cran()
  td <- .make_cytoset_test_data()

  tdr <- RunTDR(td$cs,
                .sample.var = "Sample",
                .verbose = FALSE,
                .seed = 42)

  expect_error(tdr@config$source.env$cs <- "overwrite")
})

# ======================================================================
# .get_sample_matrix correctness
# ======================================================================

test_that(".get_sample_matrix works for cytoset backend", {
  skip_on_cran()
  td <- .make_cytoset_test_data(n_cells = 50, n_markers = 5, n_samples = 4)

  tdr <- RunTDR(td$cs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .verbose = FALSE,
                .seed = 42)

  mat <- tinydenseR:::.get_sample_matrix(NULL, tdr, 1)
  expect_equal(ncol(mat), 5L)
  expect_equal(colnames(mat), paste0("Marker", 1:5))
  expect_equal(nrow(mat), 50L)
})

# ======================================================================
# .celltype.vec support
# ======================================================================

test_that("RunTDR.cytoset supports .celltype.vec", {
  skip_on_cran()
  td <- .make_cytoset_test_data(n_cells = 100, n_markers = 5, n_samples = 4)

  n_cells <- 100L
  n_samples <- 4L
  n_total <- n_cells * n_samples
  ct_vec <- rep(c("typeA", "typeB"), length.out = n_total)
  names(ct_vec) <- rep(paste0("event_", seq_len(n_cells)), times = n_samples)

  tdr <- RunTDR(td$cs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .celltype.vec = ct_vec,
                .verbose = FALSE,
                .seed = 42)

  expect_true(length(tdr@landmark.annot$celltyping) > 0)
})

# ======================================================================
# .min.cells.per.sample filtering
# ======================================================================

test_that("RunTDR.cytoset filters tiny samples", {
  skip_on_cran()
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  # Create 3 normal samples + 1 tiny sample (5 cells)
  sample_names <- paste0("sample", 1:4)
  tmp_dir <- tempfile("cytoset_tiny_")
  dir.create(tmp_dir)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  fcs_paths <- vapply(seq_along(sample_names), function(i) {
    nc <- if (i == 4) 5L else 100L
    mat <- matrix(runif(nc * 5L), nrow = nc, ncol = 5L,
                  dimnames = list(paste0("event_", seq_len(nc)),
                                  paste0("Marker", 1:5)))
    ff <- flowCore::flowFrame(exprs = mat)
    fpath <- file.path(tmp_dir, paste0(sample_names[i], ".fcs"))
    flowCore::write.FCS(ff, filename = fpath)
    fpath
  }, character(1))

  cs <- flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(fcs_paths, sample_names))

  pd <- flowWorkspace::pData(cs)
  pd$Sample <- sample_names
  pd$group <- c("A", "B", "A", "B")
  flowWorkspace::pData(cs) <- pd

  tdr <- RunTDR(cs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .min.cells.per.sample = 10,
                .verbose = FALSE,
                .seed = 42)

  expect_false("sample4" %in% names(tdr@cells))
  expect_equal(length(tdr@cells), 3L)
})
