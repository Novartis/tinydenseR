#####
# Copyright 2025 Novartis Biomedical Research Inc.
#
# SPDX-License-Identifier: MIT
#####

# ======================================================================
# Regression tests: RunTDR vs explicit stepwise workflow parity
#
# These tests ensure that RunTDR (the convenience wrapper) produces
# numerically identical results to the explicit step-by-step workflow:
#   get.meta → get.cells → setup.tdr.obj → get.landmarks → get.graph →
#   get.map → get.embedding
#
# Root cause of historical divergence: RunTDR built .cells from table()
# which alphabetically sorts sample names, while the explicit workflow
# preserves get.meta()'s first-occurrence order.  Different sample
# orderings cause different random number sequences in landmark
# selection, propagating to all downstream quantities.
# ======================================================================

library(testthat)
library(tinydenseR)

# ── helpers ───────────────────────────────────────────────────────────

#' Build a small SingleCellExperiment whose samples are NOT in
#' alphabetical order in colData.  This is the minimal reproducer
#' for the historical sample-ordering bug.
.make_test_sce <- function(n_genes = 80, n_cells_per_sample = 60,
                           seed = 42) {
  set.seed(seed)

  # deliberately NON-alphabetical sample arrangement
  sample_labels <- c("B_R2", "A_R1", "B_R1", "A_R2")
  conditions    <- c("B",    "A",    "B",    "A")
  replicates    <- c("R2",   "R1",   "R1",   "R2")

  counts_list <- lapply(seq_along(sample_labels), function(i) {
    # simple Poisson counts
    mat <- matrix(
      rpois(n_genes * n_cells_per_sample,
            lambda = sample(5:20, n_genes, replace = TRUE)),
      nrow = n_genes,
      ncol = n_cells_per_sample,
      dimnames = list(
        paste0("Gene", seq_len(n_genes)),
        paste0(sample_labels[i], "_cell", seq_len(n_cells_per_sample))
      )
    )
    methods::as(mat, "dgCMatrix")
  })

  combined <- do.call(cbind, counts_list)

  col_data <- data.frame(
    Sample    = rep(sample_labels, each = n_cells_per_sample),
    Condition = rep(conditions, each = n_cells_per_sample),
    Replicate = rep(replicates, each = n_cells_per_sample),
    row.names = colnames(combined),
    stringsAsFactors = FALSE
  )

  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = combined),
    colData = S4Vectors::DataFrame(col_data)
  )
}

# ── SCE parity ────────────────────────────────────────────────────────

test_that("RunTDR.SCE produces identical landmarks, PCA, graph, map,
           and sample embedding as the explicit workflow", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  sce <- .make_test_sce()

  # ---- explicit workflow ----
  .meta  <- get.meta(.obj = sce, .sample.var = "Sample", .verbose = FALSE)
  .cells <- get.cells(.exprs = sce, .meta = .meta,
                      .sample.var = "Sample",
                      .verbose = FALSE)[rownames(.meta)]

  set.seed(123)
  tdr_explicit <- setup.tdr.obj(
    .cells     = .cells,
    .meta      = .meta,
    .assay.type = "RNA",
    .verbose   = FALSE
  ) |>
    get.landmarks(.nHVG = 50, .nPC = 5, .verbose = FALSE, .seed = 123) |>
    get.graph(.cl.resolution.parameter = 1, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123) |>
    get.embedding(.verbose = FALSE)

  # ---- RunTDR wrapper ----
  set.seed(123)
  sce_run <- RunTDR(
    x           = sce,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG       = 50,
    .nPC        = 5,
    .cl.resolution.parameter = 1,
    .verbose    = FALSE,
    .seed       = 123
  )
  tdr_wrapper <- GetTDR(sce_run)

  # ---- assertions (strict equality) ----

  # Sample ordering preserved

  expect_identical(
    names(tdr_explicit@cells),
    names(tdr_wrapper@cells),
    info = "sample ordering must match between explicit and RunTDR"
  )

  # Metadata rows in same order
  expect_identical(
    rownames(tdr_explicit@metadata),
    rownames(tdr_wrapper@metadata),
    info = "metadata row order must match"
  )

  # Landmark raw assay
  expect_equal(
    tdr_explicit@assay$raw,
    tdr_wrapper@assay$raw,
    tolerance = 0,
    info = "raw landmarks must be identical"
  )

  # PCA coordinates
  expect_equal(
    tdr_explicit@landmark.embed$pca$coord,
    tdr_wrapper@landmark.embed$pca$coord,
    tolerance = 0,
    info = "PCA coordinates must be identical"
  )

  # PCA rotation (loadings)
  expect_equal(
    tdr_explicit@landmark.embed$pca$rotation,
    tdr_wrapper@landmark.embed$pca$rotation,
    tolerance = 0,
    info = "PCA rotation must be identical"
  )

  # PCA HVG names
  expect_identical(
    tdr_explicit@landmark.embed$pca$HVG,
    tdr_wrapper@landmark.embed$pca$HVG,
    info = "HVG names must be identical"
  )

  # Clustering
  expect_identical(
    tdr_explicit@landmark.annot$clustering$ids,
    tdr_wrapper@landmark.annot$clustering$ids,
    info = "cluster IDs must be identical"
  )

  # SNN graph
  expect_equal(
    tdr_explicit@graphs$snn,
    tdr_wrapper@graphs$snn,
    tolerance = 0,
    info = "SNN graph must be identical"
  )

  # Density matrix (fdens)
  expect_equal(
    tdr_explicit@density$fdens,
    tdr_wrapper@density$fdens,
    tolerance = 0,
    info = "fuzzy density matrix must be identical"
  )

  # Sample embedding (PCA)
  expect_equal(
    tdr_explicit@sample.embed$pca,
    tdr_wrapper@sample.embed$pca,
    tolerance = 0,
    info = "sample PCA embedding must be identical"
  )
})

# ── Seurat parity ─────────────────────────────────────────────────────

test_that("RunTDR.Seurat produces identical results as the explicit
           workflow", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  sce <- .make_test_sce()

  counts_mat <- SummarizedExperiment::assay(sce, "counts")
  col_data   <- as.data.frame(SummarizedExperiment::colData(sce))

  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts   = counts_mat,
    meta.data = col_data
  )

  # ---- explicit workflow (from Seurat) ----
  .meta  <- get.meta(.obj = seurat_obj, .sample.var = "Sample",
                     .verbose = FALSE)
  .cells <- get.cells(.exprs = seurat_obj, .meta = .meta,
                      .sample.var = "Sample",
                      .verbose = FALSE)[rownames(.meta)]

  set.seed(123)
  tdr_explicit <- setup.tdr.obj(
    .cells     = .cells,
    .meta      = .meta,
    .assay.type = "RNA",
    .verbose   = FALSE
  ) |>
    get.landmarks(.nHVG = 50, .nPC = 5, .verbose = FALSE, .seed = 123) |>
    get.graph(.cl.resolution.parameter = 1, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123) |>
    get.embedding(.verbose = FALSE)

  # ---- RunTDR wrapper ----
  set.seed(123)
  seurat_run <- RunTDR(
    x           = seurat_obj,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG       = 50,
    .nPC        = 5,
    .cl.resolution.parameter = 1,
    .verbose    = FALSE,
    .seed       = 123
  )
  tdr_wrapper <- GetTDR(seurat_run)

  # ---- assertions ----

  expect_identical(
    names(tdr_explicit@cells),
    names(tdr_wrapper@cells),
    info = "Seurat: sample ordering must match"
  )

  expect_equal(
    tdr_explicit@assay$raw,
    tdr_wrapper@assay$raw,
    tolerance = 0,
    info = "Seurat: raw landmarks must be identical"
  )

  expect_equal(
    tdr_explicit@landmark.embed$pca$coord,
    tdr_wrapper@landmark.embed$pca$coord,
    tolerance = 0,
    info = "Seurat: PCA coordinates must be identical"
  )

  expect_identical(
    tdr_explicit@landmark.annot$clustering$ids,
    tdr_wrapper@landmark.annot$clustering$ids,
    info = "Seurat: cluster IDs must be identical"
  )

  expect_equal(
    tdr_explicit@density$fdens,
    tdr_wrapper@density$fdens,
    tolerance = 0,
    info = "Seurat: fuzzy density matrix must be identical"
  )

  expect_equal(
    tdr_explicit@sample.embed$pca,
    tdr_wrapper@sample.embed$pca,
    tolerance = 0,
    info = "Seurat: sample PCA embedding must be identical"
  )
})

# ── dgCMatrix parity ─────────────────────────────────────────────────

test_that("RunTDR.dgCMatrix produces identical results as the explicit
           workflow", {
  skip_on_cran()

  sce <- .make_test_sce()

  counts_mat <- SummarizedExperiment::assay(sce, "counts")
  cell_meta  <- as.data.frame(SummarizedExperiment::colData(sce))

  # ---- explicit workflow (from list of matrices) ----
  sample_labels <- unique(cell_meta$Sample)
  count_list <- lapply(
    stats::setNames(sample_labels, sample_labels),
    function(s) counts_mat[, cell_meta$Sample == s, drop = FALSE]
  )

  .meta_explicit <- cell_meta[, c("Sample", "Condition", "Replicate")] |>
    dplyr::distinct() |>
    as.data.frame()
  rownames(.meta_explicit) <- .meta_explicit$Sample
  # match order of count_list
  .meta_explicit <- .meta_explicit[sample_labels, , drop = FALSE]

  .cells <- get.cells.list.mat(
    .count.mat.list = count_list[rownames(.meta_explicit)],
    .meta = .meta_explicit,
    .verbose = FALSE
  )

  set.seed(123)
  tdr_explicit <- setup.tdr.obj(
    .cells     = .cells,
    .meta      = .meta_explicit,
    .assay.type = "RNA",
    .verbose   = FALSE
  ) |>
    get.landmarks(.nHVG = 50, .nPC = 5, .verbose = FALSE, .seed = 123) |>
    get.graph(.cl.resolution.parameter = 1, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123) |>
    get.embedding(.verbose = FALSE)

  # ---- RunTDR.dgCMatrix ----
  set.seed(123)
  tdr_wrapper <- RunTDR(
    x           = counts_mat,
    .cell.meta  = cell_meta,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG       = 50,
    .nPC        = 5,
    .cl.resolution.parameter = 1,
    .verbose    = FALSE,
    .seed       = 123
  )

  # ---- assertions ----

  expect_identical(
    names(tdr_explicit@cells),
    names(tdr_wrapper@cells),
    info = "dgCMatrix: sample ordering must match"
  )

  expect_equal(
    tdr_explicit@assay$raw,
    tdr_wrapper@assay$raw,
    tolerance = 0,
    info = "dgCMatrix: raw landmarks must be identical"
  )

  expect_equal(
    tdr_explicit@landmark.embed$pca$coord,
    tdr_wrapper@landmark.embed$pca$coord,
    tolerance = 0,
    info = "dgCMatrix: PCA coordinates must be identical"
  )

  expect_identical(
    tdr_explicit@landmark.annot$clustering$ids,
    tdr_wrapper@landmark.annot$clustering$ids,
    info = "dgCMatrix: cluster IDs must be identical"
  )

  expect_equal(
    tdr_explicit@density$fdens,
    tdr_wrapper@density$fdens,
    tolerance = 0,
    info = "dgCMatrix: fuzzy density matrix must be identical"
  )
})

# ── Sample ordering regression ────────────────────────────────────────

test_that("RunTDR preserves get.meta() sample ordering even when
           alphabetical order differs", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")

  sce <- .make_test_sce()

  # get.meta() preserves first-occurrence order from colData
  .meta <- get.meta(.obj = sce, .sample.var = "Sample", .verbose = FALSE)
  meta_order <- rownames(.meta)

  # Verify our test SCE has non-alphabetical colData ordering
  expect_false(
    identical(meta_order, sort(meta_order)),
    info = "test SCE must have non-alphabetical sample ordering"
  )

  # RunTDR must match this exact ordering
  sce_run <- RunTDR(
    x = sce,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG = 50,
    .nPC = 5,
    .cl.resolution.parameter = 1,
    .verbose = FALSE,
    .seed = 123
  )
  tdr <- GetTDR(sce_run)

  # Also check valid-filtered ordering: some samples might be dropped,
  # but the relative order must follow .meta, not alphabetical.
  tdr_order <- names(tdr@cells)

  # All RunTDR sample names must exist in meta_order
  expect_true(all(tdr_order %in% meta_order))

  # The relative ordering must match
  meta_positions <- match(tdr_order, meta_order)
  expect_identical(
    meta_positions,
    sort(meta_positions),
    info = "RunTDR sample ordering must be a subsequence of get.meta() ordering"
  )
})

test_that("RunTDR.Seurat preserves get.meta() sample ordering", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")

  sce <- .make_test_sce()
  counts_mat <- SummarizedExperiment::assay(sce, "counts")
  col_data   <- as.data.frame(SummarizedExperiment::colData(sce))

  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = counts_mat, meta.data = col_data
  )

  .meta <- get.meta(.obj = seurat_obj, .sample.var = "Sample",
                    .verbose = FALSE)
  meta_order <- rownames(.meta)

  seurat_run <- RunTDR(
    x = seurat_obj,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG = 50,
    .nPC  = 5,
    .cl.resolution.parameter = 1,
    .verbose = FALSE,
    .seed = 123
  )
  tdr <- GetTDR(seurat_run)

  tdr_order <- names(tdr@cells)
  meta_positions <- match(tdr_order, meta_order)
  expect_identical(
    meta_positions,
    sort(meta_positions),
    info = "Seurat RunTDR sample ordering must follow get.meta() order"
  )
})

# ── DelayedMatrix parity (if BPCells available) ──────────────────────

test_that("RunTDR.DelayedMatrix preserves get.meta() sample ordering", {
  skip_on_cran()
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("BPCells")

  sce <- .make_test_sce()
  counts_mat <- SummarizedExperiment::assay(sce, "counts")

  delayed_mat <- DelayedArray::DelayedArray(counts_mat)
  cell_meta   <- as.data.frame(SummarizedExperiment::colData(sce))

  # Get first-occurrence order
  meta_order <- cell_meta[, c("Sample", "Condition", "Replicate")] |>
    dplyr::distinct() |>
    (\(x) x$Sample)()

  tdr <- RunTDR(
    x          = delayed_mat,
    .cell.meta = cell_meta,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG      = 50,
    .nPC       = 5,
    .cl.resolution.parameter = 1,
    .verbose   = FALSE,
    .seed      = 123
  )

  tdr_order <- names(tdr@cells)
  meta_positions <- match(tdr_order, meta_order)
  expect_identical(
    meta_positions,
    sort(meta_positions),
    info = "DelayedMatrix RunTDR sample ordering must follow first-occurrence order"
  )
})

# ── SCE with DelayedMatrix assay parity ──────────────────────────────

test_that("RunTDR.SCE with DelayedMatrix assay preserves sample ordering", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("BPCells")

  sce <- .make_test_sce()

  # Replace counts with a DelayedMatrix
  counts_delayed <- DelayedArray::DelayedArray(
    SummarizedExperiment::assay(sce, "counts")
  )
  SummarizedExperiment::assay(sce, "counts") <- counts_delayed

  .meta <- get.meta(.obj = sce, .sample.var = "Sample", .verbose = FALSE)
  meta_order <- rownames(.meta)

  sce_run <- RunTDR(
    x = sce,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG = 50,
    .nPC = 5,
    .cl.resolution.parameter = 1,
    .verbose = FALSE,
    .seed = 123
  )
  tdr <- GetTDR(sce_run)

  tdr_order <- names(tdr@cells)
  meta_positions <- match(tdr_order, meta_order)
  expect_identical(
    meta_positions,
    sort(meta_positions),
    info = "Delayed-backed SCE RunTDR sample ordering must follow get.meta() order"
  )
})

# ── Package data reproducer parity ────────────────────────────────────

test_that("RunTDR on sim_trajectory_tdr SCE matches explicit workflow", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")

  data(sim_trajectory_tdr, package = "tinydenseR")
  sim <- sim_trajectory_tdr
  rm(sim_trajectory_tdr)
  sim_trajectory <- sim$sce

  # ---- explicit workflow (mirrors inst/analysis script) ----
  .meta  <- get.meta(.obj = sim_trajectory, .sample.var = "Sample",
                     .verbose = FALSE)
  .cells <- get.cells(.exprs = sim_trajectory, .meta = .meta,
                      .sample.var = "Sample",
                      .verbose = FALSE)[rownames(.meta)]

  set.seed(123)
  tdr_explicit <- setup.tdr.obj(
    .cells     = .cells,
    .meta      = .meta,
    .assay.type = "RNA",
    .verbose   = FALSE
  ) |>
    get.landmarks(.nHVG = 500, .verbose = FALSE, .seed = 123) |>
    get.graph(.cl.resolution.parameter = 10, .verbose = FALSE, .seed = 123) |>
    get.map(.verbose = FALSE, .seed = 123) |>
    get.embedding(.verbose = FALSE)

  # ---- RunTDR wrapper ----
  set.seed(123)
  sce_run <- RunTDR(
    x           = sim_trajectory,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG       = 500,
    .cl.resolution.parameter = 10,
    .verbose    = FALSE,
    .seed       = 123
  )
  tdr_wrapper <- GetTDR(sce_run)

  # ---- assertions ----
  expect_identical(
    names(tdr_explicit@cells),
    names(tdr_wrapper@cells),
    info = "sim_trajectory: sample ordering must match"
  )

  expect_equal(
    tdr_explicit@assay$raw,
    tdr_wrapper@assay$raw,
    tolerance = 0,
    info = "sim_trajectory: raw landmarks must be identical"
  )

  expect_equal(
    tdr_explicit@landmark.embed$pca$coord,
    tdr_wrapper@landmark.embed$pca$coord,
    tolerance = 0,
    info = "sim_trajectory: PCA coordinates must be identical"
  )

  expect_identical(
    tdr_explicit@landmark.annot$clustering$ids,
    tdr_wrapper@landmark.annot$clustering$ids,
    info = "sim_trajectory: cluster IDs must be identical"
  )

  expect_equal(
    tdr_explicit@density$fdens,
    tdr_wrapper@density$fdens,
    tolerance = 0,
    info = "sim_trajectory: fuzzy density matrix must be identical"
  )

  expect_equal(
    tdr_explicit@sample.embed$pca,
    tdr_wrapper@sample.embed$pca,
    tolerance = 0,
    info = "sim_trajectory: sample PCA embedding must be identical"
  )
})
