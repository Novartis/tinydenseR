#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ======================================================================
# Helper: create a flowSet + pdata for testing
# ======================================================================

.make_flowset_test_data <- function(n_cells = 100, n_markers = 5,
                                    n_samples = 4) {
  skip_if_not_installed("flowCore")

  sample_names <- paste0("sample", seq_len(n_samples))

  flowframes <- lapply(sample_names, function(sn) {
    mat <- matrix(runif(n_cells * n_markers),
                  nrow = n_cells, ncol = n_markers,
                  dimnames = list(paste0("event_", seq_len(n_cells)),
                                  paste0("Marker", seq_len(n_markers))))
    flowCore::flowFrame(exprs = mat)
  })
  names(flowframes) <- sample_names

  fs <- methods::as(flowframes, "flowSet")
  flowCore::sampleNames(fs) <- sample_names

  pd <- flowCore::pData(fs)
  pd$Sample <- sample_names
  pd$group <- rep(c("A", "B"), length.out = n_samples)
  flowCore::pData(fs) <- pd

  list(fs = fs, pdata = pd)
}

# ======================================================================
# S3 dispatch
# ======================================================================

test_that("RunTDR.flowSet S3 method exists", {
  skip_if_not_installed("flowCore")
  expect_true(is.function(getS3method("RunTDR", "flowSet")))
})

test_that("RunTDR dispatches to RunTDR.flowSet for flowSet input", {
  skip_if_not_installed("flowCore")
  td <- .make_flowset_test_data(n_cells = 50, n_markers = 5, n_samples = 4)
  # Verify class is flowSet (not cytoset)
  expect_true(inherits(td$fs, "flowSet"))
  expect_false(inherits(td$fs, "cytoset"))
})

# ======================================================================
# Input validation
# ======================================================================

test_that("RunTDR.flowSet errors when .sample.var not in pData", {
  td <- .make_flowset_test_data()
  expect_error(
    RunTDR(td$fs, .sample.var = "nonexistent"),
    "not found in pData"
  )
})

test_that("RunTDR.flowSet errors when .assay.type is not cyto", {
  td <- .make_flowset_test_data()
  expect_error(
    RunTDR(td$fs, .sample.var = "Sample", .assay.type = "RNA"),
    "must be 'cyto'"
  )
})

test_that("RunTDR.flowSet errors when .markers not in channels", {
  td <- .make_flowset_test_data()
  expect_error(
    RunTDR(td$fs, .sample.var = "Sample", .markers = c("FakeMarker")),
    "not found in channels"
  )
})

# ======================================================================
# Full pipeline
# ======================================================================

test_that("RunTDR.flowSet full pipeline returns TDRObj", {
  skip_on_cran()
  td <- .make_flowset_test_data()

  tdr <- RunTDR(td$fs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .verbose = FALSE,
                .seed = 42)

  expect_true(is.TDRObj(tdr))
  expect_equal(tdr@config$backend, "cytoset")
  expect_true(length(tdr@density$fdens) > 0)
  expect_equal(nrow(tdr@metadata), length(tdr@cells))
})

# ======================================================================
# Locked environment
# ======================================================================

test_that("RunTDR.flowSet locks source env", {
  skip_on_cran()
  td <- .make_flowset_test_data()

  tdr <- RunTDR(td$fs,
                .sample.var = "Sample",
                .verbose = FALSE,
                .seed = 42)

  expect_error(tdr@config$source.env$cs <- "overwrite")
})

# ======================================================================
# .get_sample_matrix correctness
# ======================================================================

test_that(".get_sample_matrix works for flowSet backend", {
  skip_on_cran()
  td <- .make_flowset_test_data(n_cells = 50, n_markers = 5, n_samples = 4)

  tdr <- RunTDR(td$fs,
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
# Marker subsetting
# ======================================================================

test_that("RunTDR.flowSet respects .markers subset", {
  skip_on_cran()
  td <- .make_flowset_test_data(n_cells = 100, n_markers = 5, n_samples = 4)

  tdr <- RunTDR(td$fs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:3),
                .verbose = FALSE,
                .seed = 42)

  expect_equal(tdr@config$markers, paste0("Marker", 1:3))
  mat <- tinydenseR:::.get_sample_matrix(NULL, tdr, 1)
  expect_true(all(paste0("Marker", 1:3) %in% colnames(mat)))
})

# ======================================================================
# .markers = NULL defaults to all channels
# ======================================================================

test_that("RunTDR.flowSet defaults .markers to all channels", {
  skip_on_cran()
  td <- .make_flowset_test_data(n_cells = 100, n_markers = 5, n_samples = 4)

  tdr <- RunTDR(td$fs,
                .sample.var = "Sample",
                .markers = NULL,
                .verbose = FALSE,
                .seed = 42)

  expect_equal(length(tdr@config$markers), 5L)
})

# ======================================================================
# .celltype.vec support
# ======================================================================

test_that("RunTDR.flowSet supports .celltype.vec", {
  skip_on_cran()
  td <- .make_flowset_test_data(n_cells = 100, n_markers = 5, n_samples = 4)

  n_cells <- 100L
  n_samples <- 4L
  n_total <- n_cells * n_samples
  ct_vec <- rep(c("typeA", "typeB"), length.out = n_total)
  names(ct_vec) <- rep(paste0("event_", seq_len(n_cells)), times = n_samples)

  tdr <- RunTDR(td$fs,
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

test_that("RunTDR.flowSet filters tiny samples", {
  skip_on_cran()
  skip_if_not_installed("flowCore")

  # Create 3 normal samples + 1 tiny sample (5 cells)
  sample_names <- paste0("sample", 1:4)

  flowframes <- lapply(seq_along(sample_names), function(i) {
    nc <- if (i == 4) 5L else 100L
    mat <- matrix(runif(nc * 5L), nrow = nc, ncol = 5L,
                  dimnames = list(paste0("event_", seq_len(nc)),
                                  paste0("Marker", 1:5)))
    flowCore::flowFrame(exprs = mat)
  })
  names(flowframes) <- sample_names

  fs <- methods::as(flowframes, "flowSet")
  flowCore::sampleNames(fs) <- sample_names

  pd <- flowCore::pData(fs)
  pd$Sample <- sample_names
  pd$group <- c("A", "B", "A", "B")
  flowCore::pData(fs) <- pd

  tdr <- RunTDR(fs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .min.cells.per.sample = 10,
                .verbose = FALSE,
                .seed = 42)

  expect_false("sample4" %in% names(tdr@cells))
  expect_equal(length(tdr@cells), 3L)
})

# ======================================================================
# pData metadata propagation
# ======================================================================

test_that("RunTDR.flowSet propagates pData to TDRObj metadata", {
  skip_on_cran()
  td <- .make_flowset_test_data()

  tdr <- RunTDR(td$fs,
                .sample.var = "Sample",
                .markers = paste0("Marker", 1:5),
                .verbose = FALSE,
                .seed = 42)

  expect_true("group" %in% colnames(tdr@metadata))
  expect_equal(sort(rownames(tdr@metadata)), sort(td$pdata$Sample))
})

# ======================================================================
# Backend equivalence: flowSet vs cytoset produce identical results
# ======================================================================

test_that("flowSet and cytoset produce numerically equivalent TDRObj", {
  skip_on_cran()
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspace")

  # Create shared data
  n_cells <- 80L
  n_markers <- 5L
  n_samples <- 4L
  sample_names <- paste0("sample", seq_len(n_samples))
  tmp_dir <- tempfile("equiv_test_")
  dir.create(tmp_dir)
  withr::defer(unlink(tmp_dir, recursive = TRUE))

  # Write FCS files (needed for cytoset)
  fcs_paths <- vapply(sample_names, function(sn) {
    set.seed(match(sn, sample_names))
    mat <- matrix(runif(n_cells * n_markers),
                  nrow = n_cells, ncol = n_markers,
                  dimnames = list(paste0("event_", seq_len(n_cells)),
                                  paste0("Marker", seq_len(n_markers))))
    ff <- flowCore::flowFrame(exprs = mat)
    fpath <- file.path(tmp_dir, paste0(sn, ".fcs"))
    flowCore::write.FCS(ff, filename = fpath)
    fpath
  }, character(1))

  # --- Build cytoset ---
  cs <- flowWorkspace::load_cytoset_from_fcs(
    files = stats::setNames(fcs_paths, sample_names))
  # load_cytoset_from_fcs may append .fcs to names; normalise
  flowWorkspace::sampleNames(cs) <- sample_names
  pd_cs <- flowWorkspace::pData(cs)
  pd_cs$Sample <- sample_names
  pd_cs$group <- rep(c("A", "B"), length.out = n_samples)
  flowWorkspace::pData(cs) <- pd_cs

  # --- Build flowSet from same FCS files ---
  flowframes <- lapply(stats::setNames(fcs_paths, sample_names), function(fp) {
    flowCore::read.FCS(fp)
  })
  fs <- methods::as(flowframes, "flowSet")
  flowCore::sampleNames(fs) <- sample_names
  pd_fs <- flowCore::pData(fs)
  pd_fs$Sample <- sample_names
  pd_fs$group <- rep(c("A", "B"), length.out = n_samples)
  flowCore::pData(fs) <- pd_fs

  # --- Run both ---
  tdr_cs <- RunTDR(cs,
                   .sample.var = "Sample",
                   .markers = paste0("Marker", 1:5),
                   .verbose = FALSE,
                   .seed = 42)

  tdr_fs <- RunTDR(fs,
                   .sample.var = "Sample",
                   .markers = paste0("Marker", 1:5),
                   .verbose = FALSE,
                   .seed = 42)

  # --- Compare ---
  expect_equal(names(tdr_cs@cells), names(tdr_fs@cells))
  expect_equal(tdr_cs@config$markers, tdr_fs@config$markers)
  expect_equal(nrow(tdr_cs@metadata), nrow(tdr_fs@metadata))

  # Landmark expression matrices should be numerically equivalent
  expect_equal(tdr_cs@assay$expr, tdr_fs@assay$expr, tolerance = 1e-10)

  # PCA embeddings should be numerically equivalent
  expect_equal(tdr_cs@landmark.embed$pca$embed,
               tdr_fs@landmark.embed$pca$embed, tolerance = 1e-10)

  # Density vectors should be numerically equivalent
  expect_equal(tdr_cs@density$fdens, tdr_fs@density$fdens, tolerance = 1e-10)
})
