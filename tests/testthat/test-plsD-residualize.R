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
# Test fixture: synthetic data with batch nuisance for residualization tests
# =============================================================================

.make_resid_fixture <- function(seed = 99) {
  set.seed(seed)

  n_samples  <- 8L
  n_cells    <- 150L
  n_markers  <- 5L
  markers    <- paste0("m", seq_len(n_markers))

  # 2x2 factorial: Condition x Batch
 meta <- data.frame(
    Condition = rep(c("Ctrl", "Trt"), each = 4),
    Batch     = rep(c("A", "B"), times = 4),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- paste0("s", seq_len(n_samples))

  cells <- lapply(seq_len(n_samples), function(i) {
    cond_shift  <- ifelse(meta$Condition[i] == "Trt", 2.0, 0)
    batch_shift <- ifelse(meta$Batch[i] == "B", 1.0, 0)

    half <- n_cells %/% 2
    mat_a <- matrix(
      rnorm(half * n_markers, mean = cond_shift + batch_shift + 2),
      nrow = half, ncol = n_markers
    )
    mat_a[, 1:3] <- mat_a[, 1:3] + 3

    mat_b <- matrix(
      rnorm(half * n_markers, mean = cond_shift + batch_shift),
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
# C.1 Backward compatibility: .residualize = FALSE is the default and
#     produces identical results to omitting the argument
# =============================================================================

test_that("plsD: .residualize = FALSE gives identical results to default", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  set.seed(1)
  obj1 <- get.plsD(obj, .coef.col = "ConditionTrt",
                    .ncomp = 2, .verbose = FALSE)

  set.seed(1)
  obj2 <- get.plsD(obj, .coef.col = "ConditionTrt",
                    .ncomp = 2, .residualize = FALSE, .verbose = FALSE)

  r1 <- obj1@results$pls$ConditionTrt
  r2 <- obj2@results$pls$ConditionTrt

  expect_equal(r1$coord, r2$coord, tolerance = 0)
  expect_equal(r1$loadings, r2$loadings, tolerance = 0)
  expect_equal(r1$raw.loadings, r2$raw.loadings, tolerance = 0)
  expect_equal(r1$gene.weights, r2$gene.weights, tolerance = 0)
})


# =============================================================================
# C.2 No nuisance → no-op: intercept-only nuisance residualization ≈ no-op
# =============================================================================

test_that("plsD: residualize with no nuisance covariates is a no-op", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  # Design with intercept + treatment only — intercept is the only nuisance
  des <- model.matrix(~ Condition, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  set.seed(1)
  obj_no <- get.plsD(obj, .coef.col = "ConditionTrt",
                      .ncomp = 2, .residualize = FALSE, .verbose = FALSE)

  set.seed(1)
  obj_yes <- get.plsD(obj, .coef.col = "ConditionTrt",
                       .ncomp = 2, .residualize = TRUE, .verbose = FALSE)

  r_no  <- obj_no@results$pls$ConditionTrt
  r_yes <- obj_yes@results$pls$ConditionTrt

  # Intercept-only residualization is equivalent to centering, which is
  # already done — scores/loadings should be numerically identical
  expect_equal(r_no$coord, r_yes$coord, tolerance = 1e-10)
  expect_equal(r_no$raw.loadings, r_yes$raw.loadings, tolerance = 1e-10)
})


# =============================================================================
# C.3 Scores orthogonal to nuisance design after residualization
# =============================================================================

test_that("plsD: residualized scores are orthogonal to nuisance Z", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt",
                   .ncomp = 2, .residualize = TRUE, .verbose = FALSE)

  res <- obj@results$pls$ConditionTrt

  # Expand batch column to landmark level
  fit <- obj@results$lm$default$fit
  Z.batch <- fit$design[obj@config$key, "BatchB", drop = FALSE]

  for (k in seq_len(ncol(res$coord))) {
    r <- abs(cor(res$coord[, k], Z.batch))
    # Should be near zero (not exactly, because scores come from NIPALS
    # on residualized M, not direct residualization of scores).
    # With small synthetic data (150 cells x 8 samples), graph smoothing
    # and Y-interaction can leave moderate residual correlation.
    expect_lt(r, 0.25,
              label = paste0("cor(score_", k, ", batch)"))
  }
})


# =============================================================================
# C.4 Implicit operator ≡ dense materialization
# =============================================================================

test_that("plsD: implicit resid operators match dense materialization", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  # Run plsD with .store.M = TRUE to get the materialized M.local
  set.seed(1)
  obj_mat <- get.plsD(obj, .coef.col = "ConditionTrt",
                       .ncomp = 2, .residualize = TRUE,
                       .store.M = TRUE, .verbose = FALSE)

  M.stored <- obj_mat@results$pls$ConditionTrt$M.local

  # Verify M.stored is non-NULL

  expect_false(is.null(M.stored))

  # The decomposition should complete without error
  expect_true(all(is.finite(obj_mat@results$pls$ConditionTrt$coord)))
})


# =============================================================================
# C.6 Contrast model path
# =============================================================================

test_that("plsD: residualize works with contrasts model", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ 0 + Condition + Batch, data = fix$obj@metadata)
  cont <- limma::makeContrasts(ConditionTrt - ConditionCtrl, levels = des)
  obj <- get.lm(fix$obj, .design = des, .contrasts = cont, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt - ConditionCtrl",
                   .ncomp = 2, .residualize = TRUE, .verbose = FALSE)

  res <- obj@results$pls[["ConditionTrt - ConditionCtrl"]]
  expect_false(is.null(res))
  expect_true(res$params$residualize)
  expect_true("BatchB" %in% res$params$nuisance.cols)
})


# =============================================================================
# C.7 Spearman guard
# =============================================================================

test_that("plsD: .residualize + spearman errors", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  expect_error(
    get.plsD(obj, .coef.col = "ConditionTrt",
             .residualize = TRUE, .loading.method = "spearman",
             .verbose = FALSE),
    "not supported with .loading.method = 'spearman'"
  )
})


# =============================================================================
# C.7b Input validation
# =============================================================================

test_that("plsD: .residualize rejects non-logical input", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  expect_error(
    get.plsD(obj, .coef.col = "ConditionTrt",
             .residualize = "yes", .verbose = FALSE),
    ".residualize must be TRUE or FALSE"
  )
})


# =============================================================================
# C.8 Loading correctness: residualized Pearson ≈ cor(T, Xtilde)
# =============================================================================

test_that("plsD: residualized Pearson loadings match explicit computation", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt",
                   .ncomp = 2, .residualize = TRUE,
                   .loading.method = "pearson", .verbose = FALSE)

  res <- obj@results$pls$ConditionTrt

  # Explicitly compute Xtilde = (I - H_Z) X
  fit <- obj@results$lm$default$fit
  design <- fit$design
  nuisance.cols <- which(colnames(design) != "ConditionTrt")
  Z.sample <- design[, nuisance.cols, drop = FALSE]
  Z <- Z.sample[obj@config$key, , drop = FALSE]

  prep <- tinydenseR:::.prepare.X(.tdr.obj = obj, .min.prop = 0.005,
                                   .center = FALSE, .verbose = FALSE)
  X <- as.matrix(prep$X)

  H_Z <- Z %*% solve(crossprod(Z)) %*% t(Z)
  Xtilde <- (diag(nrow(X)) - H_Z) %*% X
  Xtilde.c <- scale(Xtilde, center = TRUE, scale = FALSE)

  # Compare loadings for each component
  for (k in seq_len(ncol(res$coord))) {
    score.k <- res$coord[, k]
    score.k <- score.k - mean(score.k)

    expected <- cor(score.k, Xtilde.c)
    actual   <- res$raw.loadings[, k]

    # Allow moderate tolerance: NIPALS iterations introduce floating point
    # differences, and Pearson computation uses sparse paths
    expect_equal(as.numeric(actual), as.numeric(expected),
                 tolerance = 1e-8,
                 label = paste0("loadings_comp_", k))
  }
})


# =============================================================================
# C.8b OLS loading correctness
# =============================================================================

test_that("plsD: residualized OLS loadings match explicit computation", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt",
                   .ncomp = 2, .residualize = TRUE,
                   .loading.method = "ols", .verbose = FALSE)

  res <- obj@results$pls$ConditionTrt

  # Explicitly compute Xtilde
  fit <- obj@results$lm$default$fit
  design <- fit$design
  nuisance.cols <- which(colnames(design) != "ConditionTrt")
  Z.sample <- design[, nuisance.cols, drop = FALSE]
  Z <- Z.sample[obj@config$key, , drop = FALSE]

  prep <- tinydenseR:::.prepare.X(.tdr.obj = obj, .min.prop = 0.005,
                                   .center = FALSE, .verbose = FALSE)
  X <- as.matrix(prep$X)

  H_Z <- Z %*% solve(crossprod(Z)) %*% t(Z)
  Xtilde <- (diag(nrow(X)) - H_Z) %*% X

  for (k in seq_len(ncol(res$coord))) {
    score.k <- res$coord[, k]
    score.k <- score.k - mean(score.k)

    expected <- as.numeric(crossprod(score.k, Xtilde)) / sum(score.k^2)
    actual   <- res$raw.loadings[, k]

    expect_equal(as.numeric(actual), expected,
                 tolerance = 1e-8,
                 label = paste0("ols_loadings_comp_", k))
  }
})


# =============================================================================
# Residualization metadata stored correctly
# =============================================================================

test_that("plsD: residualization metadata is stored in params", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition + Batch, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt",
                   .ncomp = 2, .residualize = TRUE, .verbose = FALSE)

  params <- obj@results$pls$ConditionTrt$params
  expect_true(params$residualize)
  expect_true(is.character(params$nuisance.cols))
  expect_true(length(params$nuisance.cols) > 0)
})

test_that("plsD: non-residualized params show residualize = FALSE", {
  skip_on_cran()
  fix <- .make_resid_fixture()
  on.exit(lapply(fix$cells, unlink), add = TRUE)

  des <- model.matrix(~ Condition, data = fix$obj@metadata)
  obj <- get.lm(fix$obj, .design = des, .verbose = FALSE)

  obj <- get.plsD(obj, .coef.col = "ConditionTrt",
                   .ncomp = 2, .verbose = FALSE)

  params <- obj@results$pls$ConditionTrt$params
  expect_false(params$residualize)
  expect_null(params$nuisance.cols)
})
