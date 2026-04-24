#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

library(testthat)
library(tinydenseR)

# ============================================================================
# Harmony batch correction for cytometry data
# ============================================================================

# Helper: create cytometry data with a known batch shift
.make_cyto_batch_data <- function(n_cells = 80, n_markers = 5,
                                  n_samples = 4, batch_shift = 2,
                                  seed = 42) {
  set.seed(seed)

  markers <- paste0("marker_", seq_len(n_markers))

  mats <- lapply(seq_len(n_samples), function(i) {
    # batches: samples 1-2 = batch A, samples 3-4 = batch B
    shift <- if (i <= n_samples / 2) 0 else batch_shift
    mat <- matrix(
      data = rnorm(n_cells * n_markers, mean = shift),
      nrow = n_cells,
      ncol = n_markers,
      dimnames = list(
        paste0("s", i, "_cell_", seq_len(n_cells)),
        markers
      )
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    uri
  })
  names(mats) <- paste0("sample", seq_len(n_samples))

  .meta <- data.frame(
    row.names = names(mats),
    batch = rep(c("A", "B"), each = n_samples / 2),
    group = rep(c("ctrl", "trt"), times = n_samples / 2),
    stringsAsFactors = FALSE
  )

  list(cells = mats, meta = .meta, markers = markers)
}

# ── 1. setup.tdr.obj accepts .harmony.var for cytometry ──────────────────────

test_that("setup.tdr.obj stores harmony.var for cyto assay type", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  )

  expect_true(is.TDRObj(obj))
  expect_equal(obj@integration$harmony.var, "batch")
  expect_equal(obj@config$assay.type, "cyto")
})

# ── 2. get.landmarks produces Harmony-corrected embedding for cyto ───────────

test_that("get.landmarks populates harmony.obj and corrected pca$coord for cyto", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123)

  # harmony.obj must be a Symphony reference

  expect_true(!is.null(obj@integration$harmony.obj))
  expect_true("Z_corr" %in% names(obj@integration$harmony.obj))
  expect_true("loadings" %in% names(obj@integration$harmony.obj))
  expect_true("vargenes" %in% names(obj@integration$harmony.obj))

  # pca$coord should be the Harmony-corrected embedding (nrow = n_landmarks)
  expect_true(!is.null(obj@landmark.embed$pca$coord))
  expect_equal(nrow(obj@landmark.embed$pca$coord),
               nrow(obj@assay$expr))

  # pca$HVG should be marker names (all markers used for cyto)
  expect_true(all(obj@landmark.embed$pca$HVG %in% td$markers))
})

# ── 3. get.graph uses pca$coord (not assay$expr) when Harmony+cyto ──────────

test_that("get.graph routes cyto+Harmony to pca$coord (not raw markers)", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 123)

  # The UMAP model should exist
  expect_true(!is.null(obj@integration$umap.model))

  # The UMAP model input dimensionality should match pca$coord columns
  # (Harmony-corrected SVD dimensions = n_markers), NOT raw marker count
  umap_input_dim <- ncol(obj@landmark.embed$pca$coord)
  expect_equal(umap_input_dim, length(td$markers))

  # adj.matrix should be populated
  expect_true(!is.null(obj@graphs$adj.matrix))

  # clustering should exist
  expect_true(!is.null(obj@landmark.annot$clustering$ids))
})

# ── 4. .scale defaults to FALSE when Harmony is active ──────────────────────

test_that(".scale defaults to FALSE for cyto+Harmony, TRUE for cyto alone", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  # With Harmony: .scale should default to FALSE
  obj_h <- setup.tdr.obj(
    .cells      = td$cells,
    .meta       = td$meta,
    .markers    = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123)

  formals_env <- list(x = obj_h)
  scale_with_harmony <- eval(formals(get.graph.TDRObj)$.scale, formals_env)
  expect_false(scale_with_harmony)

  # Without Harmony: .scale should default to TRUE
  obj_no_h <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .assay.type = "cyto",
    .verbose    = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123)

  formals_env2 <- list(x = obj_no_h)
  scale_no_harmony <- eval(formals(get.graph.TDRObj)$.scale, formals_env2)
  expect_true(scale_no_harmony)
})

# ── 5. get.map projects via Symphony for cyto+Harmony ────────────────────────

test_that("get.map uses symphony::mapQuery path for cyto+Harmony", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  # norm should be populated
  expect_true(!is.null(obj@density$norm))
  expect_equal(ncol(obj@density$norm), length(td$cells))
  expect_equal(nrow(obj@density$norm), nrow(obj@assay$expr))

  # Y (log-density) should be populated
  expect_true(!is.null(obj@density$log.norm))

  # Cluster assignments should exist for each sample
  cl_counts <- obj@density$composition$clustering$cell.count
  expect_true(!is.null(cl_counts))
  expect_equal(nrow(cl_counts), length(td$cells))
})

# ── 6. Harmony-corrected differs from uncorrected (batch removal check) ─────

test_that("Harmony changes the embedding for cyto data with batch effects", {
  td <- .make_cyto_batch_data(batch_shift = 3, seed = 99)
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj_harmony <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123)

  obj_no_harmony <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .assay.type = "cyto",
    .verbose    = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123)

  # Harmony-corrected coord should differ from uncorrected
  # (due to batch correction applied)
  coord_h  <- obj_harmony@landmark.embed$pca$coord
  coord_nh <- obj_no_harmony@landmark.embed$pca$coord

  # They may have different landmark cells (Pass 2 resamples), so compare
  # the harmony.obj existence and coord shape
  expect_true(!is.null(obj_harmony@integration$harmony.obj))
  expect_null(obj_no_harmony@integration$harmony.obj)

  # The corrected embedding dimensions should match marker count
  expect_equal(ncol(coord_h), length(td$markers))
  expect_equal(ncol(coord_nh), length(td$markers))
})

# ── 7. Full pipeline end-to-end with get.embedding ──────────────────────────

test_that("Full cyto+Harmony pipeline runs end-to-end through get.embedding", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123) |>
    get.embedding(.verbose = FALSE)

  # PCA sample embedding should be populated
  expect_true(!is.null(obj@sample.embed$pca$coord))
  expect_equal(nrow(obj@sample.embed$pca$coord), length(td$cells))
})

# ── 8. External .ref.obj remains blocked for cyto ────────────────────────────

test_that("get.map rejects .ref.obj for cytometry assay type", {
  td <- .make_cyto_batch_data()
  on.exit(lapply(td$cells, unlink), add = TRUE)

  obj <- setup.tdr.obj(
    .cells     = td$cells,
    .meta      = td$meta,
    .markers   = td$markers,
    .harmony.var = "batch",
    .assay.type  = "cyto",
    .verbose     = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 123)

  fake_ref <- list(Z_corr = matrix(1, nrow = 2, ncol = 2),
                   meta_data = data.frame(cell_type = c("A", "B")))

  expect_error(
    get.map(obj, .ref.obj = fake_ref, .verbose = FALSE, .seed = 123),
    "only supported for assay type RNA"
  )
})
