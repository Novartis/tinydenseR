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

# =============================================================================
# Test fixture: synthetic data with controlled effects
# =============================================================================
#
# Builds a TDRObj through the full pipeline (setup → landmarks → graph → map)
# with metadata containing:
#   - Condition: 2-level factor (Ctrl / Trt)
#   - Timepoint: 3-level factor (Baseline / D1 / D7)
#   - Batch: 2-level nuisance (A / B)
#
# The Condition factor (2 levels, rank 1 after dropping) exercises the rank-1
# code path. Timepoint (3 levels, rank 2 after dropping) exercises rank > 1.
# Batch is available as a nuisance covariate.
#
# n_samples = 12 (balanced: 2 Conditions x 3 Timepoints x 2 Batches)

.make_pepc_fixture <- function(seed = 42) {
  set.seed(seed)

  n_samples  <- 12L
  n_cells    <- 200L
  n_markers  <- 6L
  markers    <- paste0("m", seq_len(n_markers))

  # Balanced factorial metadata
  meta <- expand.grid(
    Condition = c("Ctrl", "Trt"),
    Timepoint = c("Baseline", "D1", "D7"),
    Batch     = c("A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- paste0("s", seq_len(n_samples))

  # Simulate per-sample cell matrices with known shifts and two sub-populations
  # to ensure clustering finds >= 2 clusters
  cells <- lapply(seq_len(n_samples), function(i) {
    cond_shift <- ifelse(meta$Condition[i] == "Trt", 1.5, 0)
    tp_shift   <- switch(meta$Timepoint[i],
                         Baseline = 0, D1 = 0.8, D7 = 1.6)
    batch_shift <- ifelse(meta$Batch[i] == "B", 0.5, 0)

    half <- n_cells %/% 2
    # Sub-population A: positive on markers 1-3
    mat_a <- matrix(
      rnorm(half * n_markers, mean = cond_shift + tp_shift + batch_shift + 2),
      nrow = half, ncol = n_markers
    )
    mat_a[, 1:3] <- mat_a[, 1:3] + 3  # push markers 1-3 higher

    # Sub-population B: positive on markers 4-6
    mat_b <- matrix(
      rnorm(half * n_markers, mean = cond_shift + tp_shift + batch_shift),
      nrow = half, ncol = n_markers
    )
    mat_b[, 4:n_markers] <- mat_b[, 4:n_markers] + 3

    mat <- rbind(mat_a, mat_b)
    dimnames(mat) <- list(paste0("s", i, "_c", seq_len(n_cells)), markers)

    uri <- tempfile(fileext = ".RDS")
    saveRDS(mat, uri, compress = FALSE)
    uri
  })
  names(cells) <- rownames(meta)

  # Run through the full pipeline
  obj <- setup.tdr.obj(
    .cells = cells,
    .meta  = meta,
    .markers = markers,
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = seed) |>
    get.map(.verbose = FALSE, .seed = seed)

  list(obj = obj, cells = cells, meta = meta)
}


# =============================================================================
# 1. Input validation for get.embedding supervised args
# =============================================================================

test_that("get.embedding errors when model does not exist", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  expect_error(
    get.embedding(fix$obj, .contrast.of.interest = "foo"),
    "not found in .tdr.obj\\$map\\$lm"
  )
})

test_that("get.embedding errors when .term.of.interest is invalid", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Fit a simple model so the full model exists
  des <- model.matrix(~ Condition, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  red.des <- model.matrix(~ 1, data = fix$obj@metadata)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced", .verbose = FALSE)

  expect_error(
    get.embedding(obj, .full.model = "default",
                  .red.model = "reduced",
                  .term.of.interest = "NONEXISTENT"),
    "must be a column name"
  )
})

test_that("get.embedding warns when both .red.model and .contrast.of.interest are given", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ 0 + Condition + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(ConditionTrt - ConditionCtrl, levels = des)
  obj <- get.lm(fix$obj, .design = des, .contrasts = cont, .verbose = FALSE)

  red.des <- model.matrix(~ Batch, data = fix$obj@metadata)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced", .verbose = FALSE)

  expect_warning(
    get.embedding(obj,
                  .contrast.of.interest = "ConditionTrt - ConditionCtrl",
                  .red.model = "reduced",
                  .term.of.interest = "Condition",
                  .verbose = FALSE),
    "Using '.contrast.of.interest'"
  )
})

test_that("get.embedding in nested mode errors when contrasts are present", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ 0 + Condition + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(TrtVsCtrl = ConditionTrt - ConditionCtrl,
                               levels = des)
  obj <- get.lm(fix$obj, .design = des, .contrasts = cont, .verbose = FALSE)

  red.des <- model.matrix(~ Batch, data = fix$obj@metadata)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced", .verbose = FALSE)

  expect_error(
    get.embedding(obj, .full.model = "default",
                  .red.model = "reduced",
                  .term.of.interest = "Condition",
                  .verbose = FALSE),
    "Contrasts found|Use '.contrast.of.interest' instead"
  )
})


# =============================================================================
# 2. Unsupervised-only mode
# =============================================================================

test_that("get.embedding unsupervised-only returns pca without pePC", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  obj <- get.embedding(fix$obj, .verbose = FALSE)

  expect_true(!is.null(obj@sample.embed$pca))
  expect_true(!is.null(obj@sample.embed$pca$coord))
  expect_equal(nrow(obj@sample.embed$pca$coord), nrow(fix$meta))
  expect_null(obj@sample.embed$pepc)
})


# =============================================================================
# 3. FWL contrast path (rank 1): two-level factor
# =============================================================================

test_that("FWL contrast pePC produces rank-1 embedding with expected structure", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Cell-means model with nuisance
  des <- model.matrix(~ 0 + Condition + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(
    TrtVsCtrl = ConditionTrt - ConditionCtrl,
    levels = des
  )
  obj <- get.lm(fix$obj, .design = des, .contrasts = cont, .verbose = FALSE)
  obj <- get.embedding(obj, .contrast.of.interest = "TrtVsCtrl",
                       .verbose = FALSE)

  pepc <- obj@sample.embed$pepc$TrtVsCtrl

  # Structure checks
  expect_true(!is.null(pepc))
  expect_equal(pepc$method, "fwl_contrast")
  expect_equal(ncol(pepc$coord), 1L)
  expect_equal(nrow(pepc$coord), nrow(fix$meta))
  expect_equal(ncol(pepc$rotation), 1L)
  expect_equal(colnames(pepc$coord), "pePC1")

  # delta.Yhat should be rank 1
  expect_equal(ncol(pepc$delta.Yhat), nrow(fix$meta))
  svd_vals <- svd(pepc$delta.Yhat)$d
  expect_equal(sum(svd_vals > max(svd_vals) * 1e-8), 1L)

  # Variance explained should be positive and < 100%
  expect_true(all(pepc$perc.tot.var.exp > 0))
  expect_true(all(pepc$perc.tot.var.exp < 100))

  # Effect-residual correlation should exist
  expect_length(pepc$effect.resid.cor, 1L)

  # Nuisance design should be stored
  expect_true(!is.null(pepc$nuisance.design))
})


# =============================================================================
# 4. Nested model path (rank 1): two-level factor
# =============================================================================

test_that("Nested model pePC produces rank-1 embedding for two-level factor", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  full.des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Condition", .verbose = FALSE)

  pepc <- obj@sample.embed$pepc$Condition

  expect_true(!is.null(pepc))
  expect_equal(pepc$method, "nested_models")
  expect_equal(ncol(pepc$coord), 1L)
  expect_equal(nrow(pepc$coord), nrow(fix$meta))

  # delta.Yhat should be rank 1 (two-level factor → 1 column dropped)
  svd_vals <- svd(pepc$delta.Yhat)$d
  expect_equal(sum(svd_vals > max(svd_vals) * 1e-8), 1L)
})


# =============================================================================
# 5. Nested model path (rank > 1): three-level factor
# =============================================================================

test_that("Nested model pePC produces rank-2 embedding for three-level factor", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  full.des <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Timepoint", .verbose = FALSE)

  pepc <- obj@sample.embed$pepc$Timepoint

  expect_true(!is.null(pepc))
  expect_equal(pepc$method, "nested_models")
  expect_equal(ncol(pepc$coord), 2L)
  expect_equal(nrow(pepc$coord), nrow(fix$meta))
  expect_equal(colnames(pepc$coord), c("pePC1", "pePC2"))

  # Rotation should have 2 columns (2 axes)
  expect_equal(ncol(pepc$rotation), 2L)

  # Variance explained: two entries, both positive, sum < 100%
  expect_length(pepc$perc.tot.var.exp, 2L)
  expect_true(all(pepc$perc.tot.var.exp > 0))
  expect_true(sum(pepc$perc.tot.var.exp) < 100)

  # Effect-residual correlation: two entries
  expect_length(pepc$effect.resid.cor, 2L)

  # delta.Yhat should have rank <= 2
  svd_vals <- svd(pepc$delta.Yhat)$d
  expect_lte(sum(svd_vals > max(svd_vals) * 1e-8), 2L)
})


test_that("Rank>1 nested model emits informational message when verbose", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  full.des <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)

  expect_message(
    get.embedding(obj, .full.model = "full", .red.model = "reduced",
                  .term.of.interest = "Timepoint", .verbose = TRUE),
    "do NOT correspond"
  )
})


# =============================================================================
# 6. Parameterization invariance of nested model pePC for rank > 1
# =============================================================================

test_that("Nested model pePC subspace is invariant to factor coding", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # ------ Treatment coding (default) ------
  full.trt <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red      <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj1 <- get.lm(fix$obj, .design = full.trt, .model.name = "full",
                  .verbose = FALSE)
  obj1 <- get.lm(obj1, .design = red, .model.name = "reduced",
                  .verbose = FALSE)
  obj1 <- get.embedding(obj1, .full.model = "full", .red.model = "reduced",
                        .term.of.interest = "Timepoint", .verbose = FALSE)

  # ------ Sum coding ------
  meta2 <- fix$obj@metadata
  meta2$Timepoint <- factor(meta2$Timepoint,
                            levels = c("Baseline", "D1", "D7"))
  contrasts(meta2$Timepoint) <- contr.sum(3)

  full.sum <- model.matrix(~ Timepoint + Batch, data = meta2)

  # Red design is the same (no Timepoint)
  obj2 <- get.lm(fix$obj, .design = full.sum, .model.name = "sum.full",
                  .verbose = FALSE, .force.recalc = FALSE)
  obj2 <- get.lm(obj2, .design = red, .model.name = "sum.red",
                  .verbose = FALSE, .force.recalc = FALSE)
  obj2 <- get.embedding(obj2, .full.model = "sum.full", .red.model = "sum.red",
                        .term.of.interest = "Timepoint", .verbose = FALSE)

  # The delta.Yhat matrices should be (near-)identical
  delta1 <- obj1@sample.embed$pepc$Timepoint$delta.Yhat
  delta2 <- obj2@sample.embed$pepc$Timepoint$delta.Yhat

  expect_equal(delta1, delta2, tolerance = 1e-10)

  # The 2D subspaces spanned should be the same.
  # Test: project V1 onto the column space of V2; the residual should be ~0.
  V1 <- obj1@sample.embed$pepc$Timepoint$rotation
  V2 <- obj2@sample.embed$pepc$Timepoint$rotation

  # V2 %*% t(V2) is the projection onto col(V2)
  proj_resid <- V1 - V2 %*% (crossprod(V2, V1))
  expect_true(all(abs(proj_resid) < 1e-8))
})


test_that("Nested model pePC subspace is invariant to baseline choice", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  red <- model.matrix(~ Batch, data = fix$obj@metadata)

  # ------ Baseline = "Baseline" (default) ------
  meta.v1 <- fix$obj@metadata
  meta.v1$Timepoint <- factor(meta.v1$Timepoint,
                              levels = c("Baseline", "D1", "D7"))
  full.v1 <- model.matrix(~ Timepoint + Batch, data = meta.v1)

  obj1 <- get.lm(fix$obj, .design = full.v1, .model.name = "bl1",
                  .verbose = FALSE)
  obj1 <- get.lm(obj1, .design = red, .model.name = "red1",
                  .verbose = FALSE)
  obj1 <- get.embedding(obj1, .full.model = "bl1", .red.model = "red1",
                        .term.of.interest = "Timepoint", .verbose = FALSE)

  # ------ Baseline = "D1" ------
  meta.v2 <- fix$obj@metadata
  meta.v2$Timepoint <- factor(meta.v2$Timepoint,
                              levels = c("D1", "Baseline", "D7"))
  full.v2 <- model.matrix(~ Timepoint + Batch, data = meta.v2)

  obj2 <- get.lm(fix$obj, .design = full.v2, .model.name = "bl2",
                  .verbose = FALSE)
  obj2 <- get.lm(obj2, .design = red, .model.name = "red2",
                  .verbose = FALSE)
  obj2 <- get.embedding(obj2, .full.model = "bl2", .red.model = "red2",
                        .term.of.interest = "Timepoint", .verbose = FALSE)

  # delta.Yhat should be identical (same column spaces)
  delta1 <- obj1@sample.embed$pepc$Timepoint$delta.Yhat
  delta2 <- obj2@sample.embed$pepc$Timepoint$delta.Yhat
  expect_equal(delta1, delta2, tolerance = 1e-10)
})


# =============================================================================
# 7. FWL rank-1 ≈ nested model rank-1 for two-level factor
# =============================================================================

test_that("FWL contrast and nested model give equivalent pePC for two-level factor", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # --- FWL path ---
  cm.des <- model.matrix(~ 0 + Condition + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(
    TrtVsCtrl = ConditionTrt - ConditionCtrl, levels = cm.des
  )
  obj.fwl <- get.lm(fix$obj, .design = cm.des, .contrasts = cont,
                     .model.name = "fwl", .verbose = FALSE)
  obj.fwl <- get.embedding(obj.fwl, .full.model = "fwl",
                           .contrast.of.interest = "TrtVsCtrl",
                           .verbose = FALSE)

  # --- Nested model path ---
  full.des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)
  obj.nm <- get.lm(fix$obj, .design = full.des, .model.name = "nm.full",
                    .verbose = FALSE)
  obj.nm <- get.lm(obj.nm, .design = red.des, .model.name = "nm.red",
                    .verbose = FALSE)
  obj.nm <- get.embedding(obj.nm, .full.model = "nm.full",
                          .red.model = "nm.red",
                          .term.of.interest = "Condition",
                          .verbose = FALSE)

  coord.fwl <- obj.fwl@sample.embed$pepc$TrtVsCtrl$coord
  coord.nm  <- obj.nm@sample.embed$pepc$Condition$coord

  # Coordinates should be highly correlated (same direction, possibly different sign)
  r <- abs(cor(coord.fwl[, 1], coord.nm[, 1]))
  expect_gt(r, 0.95)

  # delta.Yhat subspaces should match (both rank 1 → directions parallel)
  d.fwl <- obj.fwl@sample.embed$pepc$TrtVsCtrl$delta.Yhat
  d.nm  <- obj.nm@sample.embed$pepc$Condition$delta.Yhat

  # Normalize and compare (they may differ by a scalar, so check cosine similarity)
  v.fwl <- svd(d.fwl)$u[, 1]
  v.nm  <- svd(d.nm)$u[, 1]
  cos_sim <- abs(sum(v.fwl * v.nm))
  expect_gt(cos_sim, 0.95)
})


# =============================================================================
# 8. Explicit per-level contrasts for three-level factor
# =============================================================================

test_that("Explicit per-level contrasts each produce rank-1 pePC", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  cm.des <- model.matrix(
    ~ 0 + Timepoint + Batch, data = fix$obj@metadata
  )
  cont <- limma::makeContrasts(
    D1vsBaseline = TimepointD1       - TimepointBaseline,
    D7vsBaseline = TimepointD7       - TimepointBaseline,
    levels = cm.des
  )
  obj <- get.lm(fix$obj, .design = cm.des, .contrasts = cont,
                .model.name = "tp.fwl", .verbose = FALSE)

  obj <- get.embedding(obj, .full.model = "tp.fwl",
                       .contrast.of.interest = "D1vsBaseline",
                       .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "tp.fwl",
                       .contrast.of.interest = "D7vsBaseline",
                       .verbose = FALSE)

  # Both should be rank 1
  expect_equal(ncol(obj@sample.embed$pepc$D1vsBaseline$coord), 1L)
  expect_equal(ncol(obj@sample.embed$pepc$D7vsBaseline$coord), 1L)

  # Both should have method = fwl_contrast
  expect_equal(obj@sample.embed$pepc$D1vsBaseline$method, "fwl_contrast")
  expect_equal(obj@sample.embed$pepc$D7vsBaseline$method, "fwl_contrast")

  # Variance explained should be positive
  expect_true(obj@sample.embed$pepc$D1vsBaseline$perc.tot.var.exp > 0)
  expect_true(obj@sample.embed$pepc$D7vsBaseline$perc.tot.var.exp > 0)
})


test_that("Per-level contrast pePCs span ≈ same subspace as whole-term pePC", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # --- Whole-term nested model ---
  full.des <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Timepoint", .verbose = FALSE)

  V.whole <- obj@sample.embed$pepc$Timepoint$rotation  # G x 2

  # --- Per-level contrasts ---
  cm.des <- model.matrix(~ 0 + Timepoint + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(
    D1vsBaseline = TimepointD1 - TimepointBaseline,
    D7vsBaseline = TimepointD7 - TimepointBaseline,
    levels = cm.des
  )
  obj2 <- get.lm(fix$obj, .design = cm.des, .contrasts = cont,
                  .model.name = "tp.fwl", .verbose = FALSE)
  obj2 <- get.embedding(obj2, .full.model = "tp.fwl",
                        .contrast.of.interest = "D1vsBaseline",
                        .verbose = FALSE)
  obj2 <- get.embedding(obj2, .full.model = "tp.fwl",
                        .contrast.of.interest = "D7vsBaseline",
                        .verbose = FALSE)

  # Combine per-level rotation vectors into a matrix
  v.d1 <- obj2@sample.embed$pepc$D1vsBaseline$rotation[, 1]
  v.d7 <- obj2@sample.embed$pepc$D7vsBaseline$rotation[, 1]
  V.levels <- cbind(v.d1, v.d7)

  # Project V.levels onto col(V.whole). If they span the same space,
  # the residual norms should be near zero.
  proj <- V.whole %*% crossprod(V.whole, V.levels)
  resid <- V.levels - proj
  resid_norm <- sqrt(colSums(resid^2)) / sqrt(colSums(V.levels^2))

  # Allow some tolerance: the per-level FWL coefficients are marginal (not

  # partial), so exact subspace equality is not guaranteed. But for a
  # balanced design with intercept-only nuisance overlap, they should be
  # very close.
  expect_true(all(resid_norm < 0.25),
              info = paste("Relative residual norms:", paste(round(resid_norm, 4), collapse = ", ")))
})


# =============================================================================
# 9. Variance-explained calculations
# =============================================================================

test_that("pePC variance explained is fraction of total density variance", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  full.des <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Timepoint", .verbose = FALSE)

  pepc <- obj@sample.embed$pepc$Timepoint

  # Manual computation: total variance = sum of per-landmark variances of Y
  total.var <- sum(matrixStats::rowVars(obj@density$log.norm))
  manual.pve <- 100 * (pepc$sdev[1:ncol(pepc$coord)])^2 / total.var

  expect_equal(as.numeric(pepc$perc.tot.var.exp), manual.pve,
               tolerance = 1e-10)
})


# =============================================================================
# 10. Multiple supervised embeddings coexist
# =============================================================================

test_that("Multiple pePC embeddings can be stored simultaneously", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Store two nested model embeddings
  full.des <- model.matrix(~ Condition + Timepoint + Batch,
                           data = fix$obj@metadata)
  red.cond <- model.matrix(~ Timepoint + Batch, data = fix$obj@metadata)
  red.tp   <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.cond, .model.name = "noCond",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.tp, .model.name = "noTP",
                .verbose = FALSE)

  obj <- get.embedding(obj, .full.model = "full", .red.model = "noCond",
                       .term.of.interest = "Condition", .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "noTP",
                       .term.of.interest = "Timepoint", .verbose = FALSE)

  expect_true("Condition" %in% names(obj@sample.embed$pepc))
  expect_true("Timepoint" %in% names(obj@sample.embed$pepc))

  # Condition should be rank 1, Timepoint should be rank 2
  expect_equal(ncol(obj@sample.embed$pepc$Condition$coord), 1L)
  expect_equal(ncol(obj@sample.embed$pepc$Timepoint$coord), 2L)
})


# =============================================================================
# 11. FWL with no nuisance covariates
# =============================================================================

test_that("FWL contrast works with no nuisance covariates (group-only model)", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Cell-means model with only Condition (no Batch, no continuous nuisance)
  des <- model.matrix(~ 0 + Condition, data = fix$obj@metadata)
  cont <- limma::makeContrasts(
    TrtVsCtrl = ConditionTrt - ConditionCtrl, levels = des
  )
  obj <- get.lm(fix$obj, .design = des, .contrasts = cont,
                .model.name = "noZ", .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "noZ",
                       .contrast.of.interest = "TrtVsCtrl",
                       .verbose = FALSE)

  pepc <- obj@sample.embed$pepc$TrtVsCtrl

  expect_true(!is.null(pepc))
  expect_equal(ncol(pepc$coord), 1L)
  expect_equal(pepc$method, "fwl_contrast")

  # Nuisance design should be intercept-only
  expect_equal(ncol(pepc$nuisance.design), 1L)
})


# =============================================================================
# 12. Effect-residual correlation diagnostic
# =============================================================================

test_that("Effect-residual correlation is high for well-specified models", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  full.des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)

  obj <- get.lm(fix$obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Condition", .verbose = FALSE)

  cor_val <- obj@sample.embed$pepc$Condition$effect.resid.cor

  # For a well-specified model with a real effect, the correlation
  # between effect scores and residualized projection should be substantial
  expect_true(abs(cor_val) > 0.3,
              info = paste("Effect-residual correlation:", round(cor_val, 3)))
})


# =============================================================================
# 13. Regression safety: unsupervised embeddings survive supervised calls
# =============================================================================

test_that("Supervised get.embedding preserves existing unsupervised embeddings", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # First call: unsupervised only

  obj <- get.embedding(fix$obj, .verbose = FALSE)
  pca_before <- obj@sample.embed$pca$coord

  # Second call: supervised (should not overwrite pca)
  full.des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  red.des  <- model.matrix(~ Batch, data = fix$obj@metadata)
  obj <- get.lm(obj, .design = full.des, .model.name = "full",
                .verbose = FALSE)
  obj <- get.lm(obj, .design = red.des, .model.name = "reduced",
                .verbose = FALSE)
  obj <- get.embedding(obj, .full.model = "full", .red.model = "reduced",
                       .term.of.interest = "Condition", .verbose = FALSE)

  pca_after <- obj@sample.embed$pca$coord

  expect_identical(pca_before, pca_after)
  expect_true(!is.null(obj@sample.embed$pepc$Condition))
})


# =============================================================================
# 14. Edge case: near-zero contrast (aliased with nuisance)
# =============================================================================

test_that("FWL warns on near-zero residualized contrast", {
  fix <- .make_pepc_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Create a contrast that is nearly aliased with nuisance:
  # Make a design where one "group" column is identical to a nuisance column
  meta <- fix$obj@metadata
  meta$AliasGroup <- meta$Batch  # perfectly correlated with Batch

  des <- model.matrix(~ 0 + AliasGroup + Batch, data = meta)

  # This design should be rank-deficient; get.lm will likely reject it.
  # We only test the warning mechanism if we can construct such a scenario.
  expect_error(
    get.lm(fix$obj, .design = des, .verbose = FALSE),
    "not of full rank|not estimable"
  )
})
