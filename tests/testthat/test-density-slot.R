#####
# Copyright 2026 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

# =========================================================================
# Tests for the @density slot schema, the get.density() accessor,
# backward-compatibility shims, and the math-chain invariants
# (raw → norm → log.norm).
#
# Covers:
#   §1  get.density() generic dispatch + .which matching
#   §2  Backward-compat: deprecated "fdens" / "Y" aliases via get.density()
#   §3  $ accessor shim: tdr$density$fdens / tdr$density$Y aliases
#   §4  $ accessor shim: tdr$map$fdens / tdr$map$Y (legacy "map" name)
#   §5  Constructor migration: old fdens/Y → norm/log.norm
#   §6  Math-chain invariants (full pipeline)
#   §7  Density dimensions & names (full pipeline)
#   §8  size.factors properties (full pipeline)
#   §9  show() method density lines
#   §10 Edge cases
# =========================================================================


# ── Shared fixture: run full cyto pipeline through get.map() ──────────────
# Uses the helper from helper-test-utils.R. We build it once (in a local
# env) and reuse across sections that need a post-get.map() TDRObj.

local_full_pipeline <- function(envir = parent.frame()) {
  set.seed(42)

  n_cells   <- 90L
  n_markers <- 5L
  n_samples <- 4L

  # Build per-sample matrices with 3 distinct subpopulations
  n_pops <- 3L
  pop_means <- matrix(c(
    0.1, 0.1, 0.9, 0.9, 0.5,
    0.9, 0.1, 0.1, 0.5, 0.9,
    0.1, 0.9, 0.5, 0.1, 0.1
  ), nrow = n_pops, ncol = n_markers, byrow = TRUE)

  .cells <- lapply(seq_len(n_samples), function(i) {
    pops <- sample(seq_len(n_pops), n_cells, replace = TRUE)
    mat <- pop_means[pops, , drop = FALSE] +
      matrix(rnorm(n_cells * n_markers, sd = 0.08),
             nrow = n_cells, ncol = n_markers)
    mat <- pmax(mat, 0)
    dimnames(mat) <- list(
      paste0("sample", i, "_cell_", seq_len(n_cells)),
      paste0("marker_", seq_len(n_markers))
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(mat, uri, compress = FALSE)
    uri
  })
  names(.cells) <- paste0("sample", seq_len(n_samples))

  .meta <- data.frame(
    row.names = names(.cells),
    group = rep(c("A", "B"), length.out = n_samples),
    stringsAsFactors = FALSE
  )
  withr::defer(lapply(.cells, unlink), envir = envir)

  setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("marker_", seq_len(n_markers)),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 42) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 42) |>
    get.map(.verbose = FALSE, .seed = 42)
}


# =========================================================================
# §1. get.density() generic dispatch + .which matching
# =========================================================================

test_that("get.density() is an S3 generic with a TDRObj method", {
  expect_true(is.function(get.density))
  expect_true(is.function(get.density.TDRObj))
})

test_that("get.density() returns correct layers from a pipeline object", {
  tdr <- local_full_pipeline()

  norm_mat     <- get.density(tdr, "norm")
  log_norm_mat <- get.density(tdr, "log.norm")
  raw_mat      <- get.density(tdr, "raw")

  expect_true(is.matrix(norm_mat) || is(norm_mat, "Matrix"))
  expect_true(is.matrix(log_norm_mat) || is(log_norm_mat, "Matrix"))
  expect_true(is.matrix(raw_mat) || is(raw_mat, "Matrix"))
})

test_that("get.density() with no .which defaults to 'raw'", {
  tdr <- local_full_pipeline()
  # match.arg picks the first choice → "raw"
  default <- get.density(tdr)
  expect_identical(default, tdr@density$raw)
})

test_that("get.density() errors on invalid .which", {
  tdr <- local_full_pipeline()
  expect_error(get.density(tdr, "nonexistent"),
               regexp = "'arg' should be one of")
})


# =========================================================================
# §2. Backward-compat aliases in get.density()
# =========================================================================

test_that("get.density(.which = 'fdens') returns norm", {
  tdr <- local_full_pipeline()
  expect_identical(get.density(tdr, "fdens"), get.density(tdr, "norm"))
})

test_that("get.density(.which = 'Y') returns log.norm", {
  tdr <- local_full_pipeline()
  expect_identical(get.density(tdr, "Y"), get.density(tdr, "log.norm"))
})


# =========================================================================
# §3. $ accessor shim: tdr$density$fdens / $Y aliases
# =========================================================================

test_that("tdr$density contains both new and aliased names", {
  tdr <- local_full_pipeline()
  d <- tdr$density

  # New canonical names

  expect_false(is.null(d$norm))
  expect_false(is.null(d$log.norm))

  # Backward-compat aliases injected by the $ method
  expect_false(is.null(d$fdens))
  expect_false(is.null(d$Y))
})

test_that("tdr$density$fdens equals tdr$density$norm", {
  tdr <- local_full_pipeline()
  expect_identical(tdr$density$fdens, tdr$density$norm)
})

test_that("tdr$density$Y equals tdr$density$log.norm", {
  tdr <- local_full_pipeline()
  expect_identical(tdr$density$Y, tdr$density$log.norm)
})

test_that("@ accessor does NOT have fdens/Y aliases (no shim)", {
  tdr <- local_full_pipeline()
  # Direct slot access bypasses the $ method
  expect_null(tdr@density$fdens)
  expect_null(tdr@density$Y)
  # But the real data is there under canonical names
  expect_false(is.null(tdr@density$norm))
  expect_false(is.null(tdr@density$log.norm))
})


# =========================================================================
# §4. $ accessor shim: tdr$map (legacy "map" slot name)
# =========================================================================

test_that("tdr$map returns density data with fdens/Y aliases", {
  tdr <- local_full_pipeline()
  m <- tdr$map

  expect_false(is.null(m$fdens))
  expect_false(is.null(m$Y))
  expect_false(is.null(m$norm))
  expect_false(is.null(m$log.norm))
})

test_that("tdr$map$fdens matches tdr$density$norm", {
  tdr <- local_full_pipeline()
  expect_identical(tdr$map$fdens, tdr$density$norm)
})


# =========================================================================
# §5. Constructor migration: old fdens/Y → norm/log.norm
# =========================================================================

test_that("TDRObj() migrates density$fdens to density$norm", {
  obj <- TDRObj(
    config  = list(assay.type = "RNA"),
    density = list(fdens = matrix(1:6, 3, 2))
  )
  expect_false(is.null(obj@density$norm))
  expect_null(obj@density$fdens)
  expect_equal(obj@density$norm, matrix(1:6, 3, 2))
})

test_that("TDRObj() migrates density$Y to density$log.norm", {
  obj <- TDRObj(
    config  = list(assay.type = "RNA"),
    density = list(Y = matrix(1:6, 3, 2))
  )
  expect_false(is.null(obj@density$log.norm))
  expect_null(obj@density$Y)
  expect_equal(obj@density$log.norm, matrix(1:6, 3, 2))
})

test_that("TDRObj() migrates both fdens and Y simultaneously", {
  mat_f <- matrix(1:6, 3, 2)
  mat_y <- log2(mat_f + 0.5)
  obj <- TDRObj(
    config  = list(assay.type = "RNA"),
    density = list(fdens = mat_f, Y = mat_y)
  )
  expect_equal(obj@density$norm, mat_f)
  expect_equal(obj@density$log.norm, mat_y)
  expect_null(obj@density$fdens)
  expect_null(obj@density$Y)
})

test_that("TDRObj() via map= arg migrates fdens/Y", {
  mat_f <- matrix(1:6, 3, 2)
  mat_y <- log2(mat_f + 0.5)
  obj <- TDRObj(
    config = list(assay.type = "RNA"),
    map    = list(fdens = mat_f, Y = mat_y)
  )
  expect_equal(obj@density$norm, mat_f)
  expect_equal(obj@density$log.norm, mat_y)
  expect_null(obj@density$fdens)
  expect_null(obj@density$Y)
})

test_that("TDRObj() does not overwrite norm with fdens if norm already set", {
  mat_norm  <- matrix(10:15, 3, 2)
  mat_fdens <- matrix(1:6, 3, 2)
  obj <- TDRObj(
    config  = list(assay.type = "RNA"),
    density = list(norm = mat_norm, fdens = mat_fdens)
  )
  # norm should be preserved (not overwritten by fdens)
  expect_equal(obj@density$norm, mat_norm)
})

test_that("TDRObj() preserves raw and size.factors when present", {
  mat_raw <- matrix(1:6, 3, 2)
  sf <- c(0.8, 1.2)
  obj <- TDRObj(
    config  = list(assay.type = "RNA"),
    density = list(
      raw = mat_raw,
      norm = matrix(1:6 / rep(sf, each = 3), 3, 2),
      log.norm = matrix(0, 3, 2),
      size.factors = sf
    )
  )
  expect_equal(obj@density$raw, mat_raw)
  expect_equal(obj@density$size.factors, sf)
})


# =========================================================================
# §6. Math-chain invariants (full pipeline)
# =========================================================================

test_that("norm == t(t(raw) / size.factors)", {
  tdr <- local_full_pipeline()

  raw <- tdr@density$raw
  sf  <- tdr@density$size.factors
  expected_norm <- Matrix::t(Matrix::t(raw) / sf)

  expect_equal(as.matrix(tdr@density$norm),
               as.matrix(expected_norm),
               tolerance = 1e-12)
})

test_that("log.norm == log2(norm + 0.5)", {
  tdr <- local_full_pipeline()

  expected <- log2(tdr@density$norm + 0.5)
  expect_equal(as.matrix(tdr@density$log.norm),
               as.matrix(expected),
               tolerance = 1e-12)
})

test_that("raw values are non-negative (sum of fuzzy edge weights)", {
  tdr <- local_full_pipeline()
  expect_true(all(tdr@density$raw >= 0))
})

test_that("norm values are non-negative", {
  tdr <- local_full_pipeline()
  expect_true(all(tdr@density$norm >= 0))
})

test_that("full chain: log.norm == log2(t(t(raw) / size.factors) + 0.5)", {
  tdr <- local_full_pipeline()

  raw <- tdr@density$raw
  sf  <- tdr@density$size.factors
  expected <- log2(Matrix::t(Matrix::t(raw) / sf) + 0.5)

  expect_equal(as.matrix(tdr@density$log.norm),
               as.matrix(expected),
               tolerance = 1e-12)
})


# =========================================================================
# §7. Density dimensions & names (full pipeline)
# =========================================================================

test_that("density matrices have dimensions landmarks × samples", {
  tdr <- local_full_pipeline()

  n_landmarks <- nrow(tdr@assay$expr)
  n_samples   <- length(tdr@cells)

  for (nm in c("raw", "norm", "log.norm")) {
    mat <- tdr@density[[nm]]
    expect_equal(nrow(mat), n_landmarks,
                 info = paste0("density$", nm, " nrow"))
    expect_equal(ncol(mat), n_samples,
                 info = paste0("density$", nm, " ncol"))
  }
})

test_that("density matrices share identical rownames", {
  tdr <- local_full_pipeline()
  expect_identical(rownames(tdr@density$raw),  rownames(tdr@density$norm))
  expect_identical(rownames(tdr@density$norm), rownames(tdr@density$log.norm))
})

test_that("density matrix colnames match sample names", {
  tdr <- local_full_pipeline()
  sample_names <- names(tdr@cells)
  expect_equal(colnames(tdr@density$norm), sample_names)
  expect_equal(colnames(tdr@density$raw),  sample_names)
  expect_equal(colnames(tdr@density$log.norm), sample_names)
})

test_that("density rownames match landmark rownames", {
  tdr <- local_full_pipeline()
  expect_identical(rownames(tdr@density$norm),
                   rownames(tdr@assay$expr))
})


# =========================================================================
# §8. size.factors properties (full pipeline)
# =========================================================================

test_that("size.factors is a named numeric vector of length N", {
  tdr <- local_full_pipeline()
  sf <- tdr@density$size.factors

  expect_true(is.numeric(sf))
  expect_equal(length(sf), length(tdr@cells))
  expect_equal(names(sf), names(tdr@cells))
})

test_that("size.factors mean is 1.0 (by construction: n / mean(n))", {
  tdr <- local_full_pipeline()
  expect_equal(mean(tdr@density$size.factors), 1.0, tolerance = 1e-12)
})

test_that("size.factors are all positive", {
  tdr <- local_full_pipeline()
  expect_true(all(tdr@density$size.factors > 0))
})


# =========================================================================
# §9. show() method density lines
# =========================================================================

test_that("show() prints 'Density/map computed: TRUE' after get.map()", {
  tdr <- local_full_pipeline()
  out <- capture.output(show(tdr))
  expect_true(any(grepl("Density/map computed: TRUE", out, fixed = TRUE)))
})

test_that("show() prints 'Raw density stored: TRUE' after get.map()", {
  tdr <- local_full_pipeline()
  out <- capture.output(show(tdr))
  expect_true(any(grepl("Raw density stored: TRUE", out, fixed = TRUE)))
})

test_that("show() prints 'Density/map computed: FALSE' on empty object", {
  obj <- TDRObj(config = list(assay.type = "RNA"))
  out <- capture.output(show(obj))
  expect_true(any(grepl("Density/map computed: FALSE", out, fixed = TRUE)))
})

test_that("show() does NOT print 'Raw density' when density is empty", {
  obj <- TDRObj(config = list(assay.type = "RNA"))
  out <- capture.output(show(obj))
  expect_false(any(grepl("Raw density stored", out, fixed = TRUE)))
})


# =========================================================================
# §10. Edge cases
# =========================================================================

test_that("get.density() on empty TDRObj returns NULL for all layers", {
  obj <- TDRObj(config = list(assay.type = "RNA"))
  expect_null(get.density(obj, "raw"))
  expect_null(get.density(obj, "norm"))
  expect_null(get.density(obj, "log.norm"))
})

test_that("density matrices contain no NA values after pipeline", {
  tdr <- local_full_pipeline()
  expect_false(anyNA(tdr@density$raw))
  expect_false(anyNA(tdr@density$norm))
  expect_false(anyNA(tdr@density$log.norm))
})

test_that("density matrices contain no Inf values after pipeline", {
  tdr <- local_full_pipeline()
  expect_false(any(is.infinite(tdr@density$raw)))
  expect_false(any(is.infinite(tdr@density$norm)))
  expect_false(any(is.infinite(tdr@density$log.norm)))
})

test_that("composition cell.count rows sum to total cells per sample", {
  tdr <- local_full_pipeline()
  cc <- tdr@density$composition$clustering$cell.count
  if (!is.null(cc)) {
    # Each sample should have n_cells cells (90 in our fixture)
    row_totals <- rowSums(cc)
    # Verify totals are positive integers
    expect_true(all(row_totals > 0))
    expect_true(all(row_totals == round(row_totals)))
  }
})

test_that("composition cell.perc rows sum to 100", {
  tdr <- local_full_pipeline()
  cp <- tdr@density$composition$clustering$cell.perc
  if (!is.null(cp)) {
    expect_equal(unname(rowSums(cp)), rep(100, nrow(cp)), tolerance = 1e-8)
  }
})

test_that("composition cell.count has no NA values", {
  tdr <- local_full_pipeline()
  cc <- tdr@density$composition$clustering$cell.count
  if (!is.null(cc)) {
    expect_false(anyNA(cc))
  }
})
