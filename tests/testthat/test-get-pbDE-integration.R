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
library(Matrix)

# =============================================================================
# Fixtures: Minimal TDRObj with realistic structure for end-to-end testing
# =============================================================================

#' Build a minimal TDRObj with "matrix" backend suitable for get.pbDE
#' @param assay.type "RNA" or "cyto"
#' @param n_genes Number of genes/markers
#' @param n_landmarks Number of landmarks
#' @param n_cells_per_sample Number of cells per sample
#' @param n_samples Number of samples (minimum 4 for design mode)
#' @param seed Random seed
build_pbDE_fixture <- function(assay.type = "RNA",
                               n_genes = 50,
                               n_landmarks = 10,
                               n_cells_per_sample = 60,
                               n_samples = 4,
                               seed = 42) {
  set.seed(seed)

  sample_names <- paste0("sample", seq_len(n_samples))
  gene_names   <- paste0("gene", seq_len(n_genes))
  n_total      <- n_cells_per_sample * n_samples

  # --- Build source expression matrix ---
  if (assay.type == "RNA") {
    # RNA: genes x cells
    src_mat <- matrix(
      rpois(n_genes * n_total, lambda = 5),
      nrow = n_genes, ncol = n_total,
      dimnames = list(gene_names, paste0("cell", seq_len(n_total)))
    )
  } else {
    # Cytometry: cells x markers
    src_mat <- matrix(
      rnorm(n_total * n_genes, mean = 3, sd = 1),
      nrow = n_total, ncol = n_genes,
      dimnames = list(paste0("cell", seq_len(n_total)), gene_names)
    )
  }

  # --- Per-sample cell indices ---
  cells_list <- stats::setNames(
    lapply(seq_len(n_samples), function(i) {
      start <- (i - 1L) * n_cells_per_sample + 1L
      end   <- i * n_cells_per_sample
      start:end
    }),
    nm = sample_names
  )

  # --- Fuzzy graphs: cells x landmarks sparse matrices with positive values ---
  fuzzy_graphs <- stats::setNames(
    lapply(seq_len(n_samples), function(i) {
      fg <- abs(rsparsematrix(
        nrow = n_cells_per_sample, ncol = n_landmarks, density = 0.4
      ))
      # Ensure no all-zero rows (every cell mapped to at least 1 landmark)
      zero_rows <- which(rowSums(fg) == 0)
      for (r in zero_rows) {
        fg[r, sample.int(n_landmarks, 1)] <- runif(1, 0.1, 1)
      }
      fg
    }),
    nm = sample_names
  )

  # --- Clustering labels: assign each landmark to one of 3 clusters ---
  cluster_ids <- factor(paste0("C", rep(1:3, length.out = n_landmarks)))

  # --- Per-cell clustering (transfer from landmarks) ---
  cell_clustering <- stats::setNames(
    lapply(seq_len(n_samples), function(i) {
      # Assign each cell to the cluster of its highest-weight landmark
      max_lm <- apply(fuzzy_graphs[[i]], 1, which.max)
      cluster_ids[max_lm]
    }),
    nm = sample_names
  )

  # --- Assay (landmark expression): L x genes ---
  assay_expr <- matrix(
    rnorm(n_landmarks * n_genes, mean = 5, sd = 2),
    nrow = n_landmarks, ncol = n_genes,
    dimnames = list(paste0("lm", seq_len(n_landmarks)), gene_names)
  )

  # --- Metadata ---
  meta <- data.frame(
    row.names = sample_names,
    n.cells   = rep(n_cells_per_sample, n_samples),
    group     = rep(c("Treatment", "Control"), length.out = n_samples),
    stringsAsFactors = FALSE
  )

  # --- Build TDRObj ---
  tdr <- TDRObj(
    cells    = cells_list,
    metadata = meta,
    config   = list(
      assay.type = assay.type,
      backend    = "matrix",
      sampling   = list(
        n.cells = stats::setNames(
          rep(n_cells_per_sample, n_samples), sample_names
        )
      ),
      source.env = list(mat = src_mat)
    ),
    assay = list(expr = assay_expr),
    landmark.annot = list(
      clustering = list(ids = cluster_ids)
    ),
    cellmap = list(
      fuzzy.graphs = fuzzy_graphs,
      clustering   = list(ids = cell_clustering)
    ),
    results = list()
  )

  return(tdr)
}

# =============================================================================
# DESIGN MODE — RNA
# =============================================================================

test_that("get.pbDE design mode runs end-to-end for RNA", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 123)

  design <- model.matrix(~ group, data = tdr@metadata)

  result <- get.pbDE(tdr, .design = design, .verbose = FALSE)

  # Results stored in correct slot

  expect_true(!is.null(result@results$pb$default$all))
  de <- result@results$pb$default$all

  # limma fit object structure

  expect_true("coefficients" %in% names(de))
  expect_true("p.value" %in% names(de))
  expect_true("adj.p" %in% names(de))
  expect_true("id.idx" %in% names(de))
  expect_true("n.pseudo" %in% names(de))

  # Dimensions: genes (filtered) x coefficients
  expect_equal(ncol(de$coefficients), ncol(design))
  expect_true(nrow(de$coefficients) > 0)

  # adj.p has same dimensions as p.value
  expect_equal(dim(de$adj.p), dim(de$p.value))
})

test_that("get.pbDE design mode produces valid statistics", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 6, seed = 456)

  design <- model.matrix(~ group, data = tdr@metadata)
  result <- get.pbDE(tdr, .design = design, .verbose = FALSE)
  de <- result@results$pb$default$all

  # p-values in [0, 1]
  expect_true(all(de$p.value >= 0 & de$p.value <= 1, na.rm = TRUE))
  expect_true(all(de$adj.p >= 0 & de$adj.p <= 1, na.rm = TRUE))

  # adj.p >= p.value (FDR correction inflates)
  expect_true(all(de$adj.p >= de$p.value - .Machine$double.eps, na.rm = TRUE))
})

test_that("get.pbDE design mode with .id subsets to population", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 789)

  design <- model.matrix(~ group, data = tdr@metadata)

  result <- get.pbDE(tdr, .design = design, .id = "C1",
                     .result.name = "C1_only", .verbose = FALSE)

  expect_true(!is.null(result@results$pb$default$C1_only))
  de <- result@results$pb$default$C1_only

  # Cell indices should be a subset
  expect_true(all(lengths(de$id.idx) > 0))
  expect_true(all(lengths(de$id.idx) < 60))  # less than all cells
})

test_that("get.pbDE design mode with contrasts", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 101)

  # Use a two-group design without intercept for clean contrast specification
  grp <- factor(tdr@metadata$group)
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  contrasts <- limma::makeContrasts(Treatment - Control, levels = design)

  result <- get.pbDE(tdr, .design = design, .contrasts = contrasts,
                     .verbose = FALSE)

  de <- result@results$pb$default$all
  expect_true(nrow(de$coefficients) > 0)
})

test_that("get.pbDE design mode .force.recalc works", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 202)
  design <- model.matrix(~ group, data = tdr@metadata)

  result <- get.pbDE(tdr, .design = design, .verbose = FALSE)

  # Running again without force should error
  expect_error(
    get.pbDE(result, .design = design, .verbose = FALSE),
    "already exist"
  )

  # With force recalc should succeed
  result2 <- get.pbDE(result, .design = design,
                      .force.recalc = TRUE, .verbose = FALSE)
  expect_true(!is.null(result2@results$pb$default$all))
})

# =============================================================================
# DESIGN MODE — Cytometry
# =============================================================================

test_that("get.pbDE design mode runs end-to-end for cytometry", {
  tdr <- build_pbDE_fixture(assay.type = "cyto", n_samples = 4, seed = 303)

  design <- model.matrix(~ group, data = tdr@metadata)

  result <- get.pbDE(tdr, .design = design, .verbose = FALSE)

  de <- result@results$pb$default$all
  expect_true("coefficients" %in% names(de))
  expect_true("adj.p" %in% names(de))
  expect_true(nrow(de$coefficients) > 0)

  # p-values valid
  expect_true(all(de$p.value >= 0 & de$p.value <= 1, na.rm = TRUE))
})

# =============================================================================
# MARKER MODE — RNA
# =============================================================================

test_that("get.pbDE marker mode runs end-to-end for RNA (vs all)", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 404)

  result <- get.pbDE(tdr, .mode = "marker", .id = "C1", .verbose = FALSE)

  # Results stored in marker slot
  expect_true(!is.null(result@results$marker))
  de <- result@results$marker$default$C1_vs_all

  expect_true("coefficients" %in% names(de))
  expect_true("adj.p" %in% names(de))
  expect_true("id1.idx" %in% names(de))
  expect_true("id2.idx" %in% names(de))
  expect_true("n.pseudo1" %in% names(de))
  expect_true("n.pseudo2" %in% names(de))

  # Valid statistics
  expect_true(all(de$p.value >= 0 & de$p.value <= 1, na.rm = TRUE))
})

test_that("get.pbDE marker mode runs with explicit .id2", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 505)

  result <- get.pbDE(tdr, .mode = "marker", .id = "C1", .id2 = "C2",
                     .verbose = FALSE)

  de <- result@results$marker$default$C1_vs_C2
  expect_true(!is.null(de))
  expect_true(nrow(de$coefficients) > 0)
})

test_that("get.pbDE marker mode with multi-cluster .id", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 606)

  result <- get.pbDE(tdr, .mode = "marker",
                     .id = c("C1", "C2"), .id2 = "C3",
                     .verbose = FALSE)

  de <- result@results$marker$default$`C1_C2_vs_C3`
  expect_true(!is.null(de))
})

# =============================================================================
# MARKER MODE — Cytometry
# =============================================================================

test_that("get.pbDE marker mode runs end-to-end for cytometry", {
  tdr <- build_pbDE_fixture(assay.type = "cyto", n_samples = 4, seed = 707)

  result <- get.pbDE(tdr, .mode = "marker", .id = "C1", .verbose = FALSE)

  de <- result@results$marker$default$C1_vs_all
  expect_true(!is.null(de))
  expect_true("coefficients" %in% names(de))
  expect_true(all(de$p.value >= 0 & de$p.value <= 1, na.rm = TRUE))
})

test_that("get.pbDE marker mode cytometry with explicit .id2", {
  tdr <- build_pbDE_fixture(assay.type = "cyto", n_samples = 4, seed = 808)

  result <- get.pbDE(tdr, .mode = "marker", .id = "C1", .id2 = "C2",
                     .verbose = FALSE)

  de <- result@results$marker$default$C1_vs_C2
  expect_true(!is.null(de))
  expect_true(nrow(de$coefficients) > 0)
})

# =============================================================================
# PROGRESS BAR (.verbose = TRUE)
# =============================================================================

test_that("get.pbDE design mode emits progress bar when .verbose = TRUE (RNA)", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 901)
  design <- model.matrix(~ group, data = tdr@metadata)

  output <- capture.output(
    result <- get.pbDE(tdr, .design = design, .verbose = TRUE),
    type = "output"
  )

  # Progress bar uses cat() which goes to stdout
  expect_true(any(grepl("\\[=+>?\\s*\\]", output)) ||
              any(grepl("100\\.0%", output)) ||
              any(grepl("done in", output)))
})

test_that("get.pbDE design mode emits stage header when .verbose = TRUE", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 902)
  design <- model.matrix(~ group, data = tdr@metadata)

  msgs <- capture.output(
    result <- get.pbDE(tdr, .design = design, .verbose = TRUE),
    type = "message"
  )

  expect_true(any(grepl("aggregating pseudobulks", msgs)))
})

test_that("get.pbDE marker mode emits two progress bars (group 1 and 2)", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 903)

  msgs <- capture.output(
    result <- get.pbDE(tdr, .mode = "marker", .id = "C1", .verbose = TRUE),
    type = "message"
  )

  expect_true(any(grepl("group 1", msgs)))
  expect_true(any(grepl("group 2", msgs)))
})

test_that("get.pbDE silent when .verbose = FALSE", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 904)
  design <- model.matrix(~ group, data = tdr@metadata)

  output <- capture.output(
    result <- get.pbDE(tdr, .design = design, .verbose = FALSE),
    type = "output"
  )
  msgs <- capture.output(
    result <- get.pbDE(tdr, .design = design, .verbose = FALSE,
                       .force.recalc = TRUE),
    type = "message"
  )

  expect_equal(length(output), 0)
  expect_equal(length(msgs), 0)
})

# =============================================================================
# RESULT NAMING
# =============================================================================

test_that("get.pbDE auto-generates .result.name in design mode", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 1001)
  design <- model.matrix(~ group, data = tdr@metadata)

  # Without .id → "all"
  result <- get.pbDE(tdr, .design = design, .verbose = FALSE)
  expect_true("all" %in% names(result@results$pb$default))

  # With .id → joined cluster names
  result2 <- get.pbDE(tdr, .design = design, .id = "C1", .verbose = FALSE)
  expect_true("C1" %in% names(result2@results$pb$default))
})

test_that("get.pbDE respects custom .model.name", {

  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 1002)
  design <- model.matrix(~ group, data = tdr@metadata)

  result <- get.pbDE(tdr, .design = design,
                     .model.name = "model_v2", .verbose = FALSE)
  expect_true("model_v2" %in% names(result@results$pb))
})

# =============================================================================
# EDGE CASES & ROBUSTNESS
# =============================================================================

test_that("get.pbDE handles sample outlier removal gracefully (RNA design)", {
  # Create fixture where one sample has very few cells in the target population
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 5, seed = 1101)

  # Artificially make sample1 have almost no cells in C1
  tdr@cellmap$clustering$ids$sample1 <- factor(
    rep("C2", length(tdr@cellmap$clustering$ids$sample1)),
    levels = c("C1", "C2", "C3")
  )

  design <- model.matrix(~ group, data = tdr@metadata)

  # Should warn about outlier removal but still complete
  expect_warning(
    result <- get.pbDE(tdr, .design = design, .id = "C1", .verbose = FALSE),
    "removed"
  )
  expect_true(!is.null(result@results$pb$default$C1))
})

test_that("get.pbDE marker mode with .id.idx works", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 1201)

  # Use landmark indices 1:4 as group 1
  result <- get.pbDE(tdr, .mode = "marker", .id.idx = 1:4, .verbose = FALSE)

  de <- result@results$marker$default[[1]]
  expect_true(!is.null(de))
  expect_true("coefficients" %in% names(de))
})

test_that("get.pbDE reproducibility: same inputs yield same results", {
  tdr <- build_pbDE_fixture(assay.type = "RNA", n_samples = 4, seed = 1301)
  design <- model.matrix(~ group, data = tdr@metadata)

  result1 <- get.pbDE(tdr, .design = design, .verbose = FALSE)
  result2 <- get.pbDE(tdr, .design = design, .verbose = FALSE,
                      .result.name = "all2")

  de1 <- result1@results$pb$default$all
  de2 <- result2@results$pb$default$all2

  expect_equal(de1$coefficients, de2$coefficients)
  expect_equal(de1$p.value, de2$p.value)
})
