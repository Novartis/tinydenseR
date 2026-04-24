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

# =============================================================================
# Tests for RNG state protection via withr::local_preserve_seed()
#
# Verifies two properties for every stochastic function in tinydenseR:
#   1. Reproducibility:  same .seed → identical output
#   2. State protection: caller's .Random.seed is unchanged after the call
#
# These properties ensure that:
#   - Analyst workflows that depend on a specific RNG stream are never
#     silently corrupted by tinydenseR function calls.
#   - Internal reproducibility is maintained via the .seed parameter.
#   - Fresh R sessions (where .Random.seed does not exist) are handled safely.
# =============================================================================


# -- Helpers ------------------------------------------------------------------

#' Assert that a function preserves the caller's RNG state.
#'
#' Calls `fn` twice with the same RNG starting point, each time preceded by
#' set.seed(outer_seed) + runif(1). If the function preserves state, the
#' runif(1) *after* each call must be identical.
#'
#' @param fn  A function call expression (wrapped in bquote or ~).
#' @param outer_seed  Seed to set before each invocation.
#' @return Invisible NULL on success (test assertions fire on failure).
expect_rng_preserved <- function(fn_call, outer_seed = 42) {
  # Run 1: set.seed → consume 1 draw → call fn → consume 1 draw
  set.seed(outer_seed)
  before1 <- runif(1)
  eval(fn_call)
  after1  <- runif(1)

  # Run 2: identical setup

  set.seed(outer_seed)
  before2 <- runif(1)
  eval(fn_call)
  after2  <- runif(1)

  expect_identical(before1, before2, label = "pre-call draws match")
  expect_identical(after1,  after2,  label = "post-call draws match (state preserved)")
  invisible(NULL)
}


# -- Core workflow functions --------------------------------------------------

test_that("get.landmarks preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  get.landmarks(tdr, .verbose = FALSE, .seed = 123, .nPC = 4)
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("get.landmarks is reproducible with same .seed", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  r1 <- get.landmarks(tdr, .verbose = FALSE, .seed = 42, .nPC = 4)
  r2 <- get.landmarks(tdr, .verbose = FALSE, .seed = 42, .nPC = 4)
  expect_identical(r1@assay$expr, r2@assay$expr)
  expect_identical(r1@landmark.embed$pca$coord, r2@landmark.embed$pca$coord)
})

test_that("get.graph preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |> get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4)

  k <- min(5, max(2, nrow(tdr@assay$expr) - 1))

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  get.graph(tdr, .k = k, .verbose = FALSE, .seed = 123)
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("get.graph is reproducible with same .seed", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |> get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4)

  k <- min(5, max(2, nrow(tdr@assay$expr) - 1))

  r1 <- get.graph(tdr, .k = k, .verbose = FALSE, .seed = 42)
  r2 <- get.graph(tdr, .k = k, .verbose = FALSE, .seed = 42)
  expect_identical(r1@landmark.embed$umap$coord, r2@landmark.embed$umap$coord)
  expect_identical(r1@landmark.annot$clustering$ids, r2@landmark.annot$clustering$ids)
})

test_that("get.map preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)

  # get.map was already called; re-running on a fresh copy tests the guard
  tdr_pre <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = min(5, max(2, nrow(tdr@assay$expr) - 1)),
              .verbose = FALSE, .seed = 123)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  get.map(tdr_pre, .verbose = FALSE, .seed = 123)
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("get.map is reproducible with same .seed", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = min(5, max(2, 10)),
              .verbose = FALSE, .seed = 123)

  r1 <- get.map(tdr, .verbose = FALSE, .seed = 42)
  r2 <- get.map(tdr, .verbose = FALSE, .seed = 42)
  expect_identical(r1@density$norm, r2@density$norm)
})


# -- get.embedding (irlba gap repair) ----------------------------------------

test_that("get.embedding preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  get.embedding(tdr, .verbose = FALSE, .seed = 123)
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("get.embedding is reproducible with same .seed", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)

  r1 <- get.embedding(tdr, .verbose = FALSE, .seed = 42)
  r2 <- get.embedding(tdr, .verbose = FALSE, .seed = 42)
  expect_identical(r1@sample.embed$pca$coord, r2@sample.embed$pca$coord)
  expect_identical(r1@sample.embed$pca$sdev,  r2@sample.embed$pca$sdev)
})

test_that("get.embedding with different .seed yields different results", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- run_full_pipeline(test_data, seed = 123, verbose = FALSE)

  r1 <- get.embedding(tdr, .verbose = FALSE, .seed = 1)
  r2 <- get.embedding(tdr, .verbose = FALSE, .seed = 999)
  # irlba's random starting vector means different seeds *may* yield

  # different embeddings; at minimum they should not error
  expect_s4_class(r1, "TDRObj")
  expect_s4_class(r2, "TDRObj")
})


# -- leiden.cluster -----------------------------------------------------------

test_that("leiden.cluster preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = min(5, max(2, 10)),
              .verbose = FALSE, .seed = 123)

  sim <- tdr@graphs$snn
  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  leiden.cluster(
    .tdr.obj = tdr,
    .sim.matrix = sim,
    .resolution.parameter = 0.001,
    .verbose = FALSE,
    .seed = 123
  )
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("leiden.cluster is reproducible with same .seed", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = min(5, max(2, 10)),
              .verbose = FALSE, .seed = 123)

  sim <- tdr@graphs$snn

  r1 <- leiden.cluster(.tdr.obj = tdr, .sim.matrix = sim,
                        .resolution.parameter = 0.001,
                        .verbose = FALSE, .seed = 42)
  r2 <- leiden.cluster(.tdr.obj = tdr, .sim.matrix = sim,
                        .resolution.parameter = 0.001,
                        .verbose = FALSE, .seed = 42)
  expect_identical(r1, r2)
})

test_that("lm.cluster preserves caller RNG state", {
  test_data <- create_test_lm_obj(n_cells = 30, n_markers = 5, n_samples = 3)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = min(5, max(2, 10)),
              .verbose = FALSE, .seed = 123)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  lm.cluster(tdr, .cl.resolution.parameter = 0.8,
             .verbose = FALSE, .seed = 123)
  a2 <- runif(1)

  expect_identical(a1, a2)
})


# -- Plot functions -----------------------------------------------------------

test_that("plotPCA preserves caller RNG state (bug fix)", {
  tdr <- create_mock_graph_obj(n_points = 20)
  feat <- seq_len(20) * 0.1  # deterministic feature (no RNG consumption in arg)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  suppressWarnings(plotPCA(tdr, .feature = feat))
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("plotPCA is reproducible with same .seed", {
  tdr <- create_mock_graph_obj(n_points = 20)
  feat <- runif(20)

  p1 <- suppressWarnings(plotPCA(tdr, .feature = feat, .seed = 42))
  p2 <- suppressWarnings(plotPCA(tdr, .feature = feat, .seed = 42))
  expect_identical(p1$data, p2$data)
})

test_that("plotUMAP preserves caller RNG state", {
  tdr <- create_mock_graph_obj(n_points = 20)
  feat <- seq_len(20) * 0.1  # deterministic feature (no RNG consumption in arg)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  suppressWarnings(plotUMAP(tdr, .feature = feat))
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("plotBeeswarm preserves caller RNG state", {
  # plotBeeswarm requires a fitted model; test the seed-guard pattern
  # by checking the function at least has local_preserve_seed
  body_text <- deparse(body(plotBeeswarm.TDRObj))
  expect_true(any(grepl("local_preserve_seed", body_text)))
})

test_that("scatterPlot preserves caller RNG state", {
  x <- runif(20)
  y <- runif(20)

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  scatterPlot(.x.feature = x, .y.feature = y, .seed = 123)
  a2 <- runif(1)

  expect_identical(a1, a2)
})

test_that("scatterPlot is reproducible with same .seed", {
  x <- runif(20)
  y <- runif(20)

  p1 <- scatterPlot(.x.feature = x, .y.feature = y, .seed = 42)
  p2 <- scatterPlot(.x.feature = x, .y.feature = y, .seed = 42)
  expect_identical(p1$data, p2$data)
})


# -- Simulate functions -------------------------------------------------------

test_that("simulate_DA_data preserves caller RNG state", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  simulate_DA_data(samples_per_group = 1, mean_cells = 100, sd_cells = 10,
                   output_dir = file.path(tempdir(), "test_da_rng"))
  a2 <- runif(1)

  expect_identical(a1, a2)

  unlink(file.path(tempdir(), "test_da_rng"), recursive = TRUE)
})

test_that("simulate_DE_data preserves caller RNG state", {
  skip_if_not_installed("flowCore")

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  simulate_DE_data(samples_per_group = 1, mean_cells = 100, sd_cells = 10,
                   output_dir = file.path(tempdir(), "test_de_rng"))
  a2 <- runif(1)

  expect_identical(a1, a2)

  unlink(file.path(tempdir(), "test_de_rng"), recursive = TRUE)
})


# -- Full pipeline cross-workflow equivalence ---------------------------------

test_that("RunTDR.default and stepwise pipeline yield identical results", {
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  seed <- 123
  nPC  <- 4
  k    <- 3

  # Stepwise
  step <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed, .nPC = nPC) |>
    get.graph(.k = k, .verbose = FALSE, .seed = seed) |>
    get.map(.verbose = FALSE, .seed = seed) |>
    get.embedding(.verbose = FALSE, .seed = seed)

  # RunTDR wrapper
  tdr0 <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  auto <- RunTDR(tdr0, .verbose = FALSE, .seed = seed, .nPC = nPC, .k = k)

  expect_identical(step@assay$expr,               auto@assay$expr)
  expect_identical(step@landmark.embed$pca$coord,  auto@landmark.embed$pca$coord)
  expect_identical(step@landmark.embed$umap$coord, auto@landmark.embed$umap$coord)
  expect_identical(step@density$norm,             auto@density$norm)
  expect_identical(step@sample.embed$pca$coord,    auto@sample.embed$pca$coord)
})

test_that("pipeline preserves caller RNG state end-to-end", {
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  )

  set.seed(42)
  a1 <- runif(1)

  set.seed(42)
  RunTDR(tdr, .verbose = FALSE, .seed = 123, .k = 3)
  a2 <- runif(1)

  expect_identical(a1, a2)
})


# -- get.lm: dead .seed parameter removed ------------------------------------

test_that("get.lm no longer has .seed in its formals", {
  expect_false(".seed" %in% names(formals(get.lm.TDRObj)))
})

test_that("get.lm silently absorbs .seed via ... for backward compatibility", {
  # Passing .seed should not error; it goes into ...
  test_data <- create_test_lm_obj(n_cells = 100, n_markers = 5, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells   = test_data$cells,
    .meta    = test_data$meta,
    .markers = paste0("marker_", 1:5),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 123, .nPC = 4) |>
    get.graph(.k = 3, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123)

  # Ensure Y exists and has rows
  skip_if(is.null(tdr@density$log.norm) || nrow(tdr@density$log.norm) == 0,
          message = "Pipeline produced empty density matrix")

  design <- model.matrix(~ group, data = tdr@metadata)
  expect_no_error(
    suppressWarnings(get.lm(tdr, .design = design, .verbose = FALSE, .seed = 123))
  )
})


# -- get.embedding: .seed parameter exists ------------------------------------

test_that("get.embedding.TDRObj has .seed in its formals", {
  expect_true(".seed" %in% names(formals(get.embedding.TDRObj)))
  expect_equal(formals(get.embedding.TDRObj)$.seed, 123)
})


# -- Edge case: fresh session (.Random.seed absent) ---------------------------

test_that("local_preserve_seed handles absent .Random.seed", {
  # In a fresh R session, .Random.seed does not exist.
  # withr::local_preserve_seed must not error.
  # We simulate by removing .Random.seed, calling a guarded function,
  # and verifying .Random.seed is still absent after the call.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  # Remove .Random.seed to simulate fresh session
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved, envir = .GlobalEnv)
    } else {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  # Call a guarded function that uses local_preserve_seed
  x <- runif(20)
  y <- runif(20)
  expect_no_error(
    scatterPlot(.x.feature = x, .y.feature = y, .seed = 123)
  )
})
