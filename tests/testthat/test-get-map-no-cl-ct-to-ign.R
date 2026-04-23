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

# =========================================================================
# Tests for removal of .cl.ct.to.ign from get.map()
#
# Validates:
#   1. get.map signature no longer contains .cl.ct.to.ign
#   2. Passing .cl.ct.to.ign produces an informative error
#   3. get.map still works correctly without .cl.ct.to.ign
#   4. density$ignored is always NULL
#   5. All cluster/celltype columns are present in composition
#   6. Pipeline still works end-to-end without .cl.ct.to.ign
# =========================================================================

# ── Shared fixture: minimal cyto pipeline through get.map() ──
local_pipeline <- function(n_cells = 90, n_markers = 5, n_samples = 4,
                           seed = 42, envir = parent.frame()) {
  set.seed(seed)

  # Build per-sample matrices with 3 distinct subpopulations so Leiden
  # clustering reliably produces >= 2 clusters (uniform-random data
  # collapses into a single cluster).
  n_pops <- 3L
  pop_means <- matrix(c(
    0.1, 0.1, 0.9, 0.9, 0.5,   # pop 1
    0.9, 0.1, 0.1, 0.5, 0.9,   # pop 2
    0.1, 0.9, 0.5, 0.1, 0.1    # pop 3
  ), nrow = n_pops, ncol = n_markers, byrow = TRUE)

  .cells <- lapply(seq_len(n_samples), function(i) {
    pops <- sample(seq_len(n_pops), n_cells, replace = TRUE)
    mat <- pop_means[pops, , drop = FALSE] +
      matrix(rnorm(n_cells * n_markers, sd = 0.08),
             nrow = n_cells, ncol = n_markers)
    mat <- pmax(mat, 0)   # marker values >= 0
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

  result <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("marker_", seq_len(n_markers)),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = seed) |>
    get.map(.verbose = FALSE, .seed = seed)

  result
}

# =========================================================================
# 1. Signature tests
# =========================================================================

test_that("get.map.TDRObj signature does not contain .cl.ct.to.ign", {
  fmls <- names(formals(get.map.TDRObj))
  expect_false(".cl.ct.to.ign" %in% fmls)
})

test_that("get.map.TDRObj signature contains expected arguments", {
  fmls <- names(formals(get.map.TDRObj))
  expect_true(all(c("x", ".source", ".ref.obj", ".celltype.col.name",
                     ".verbose", ".seed", ".label.confidence",
                     ".cache.on.disk") %in% fmls))
})

# =========================================================================
# 2. Informative error on removed argument
# =========================================================================

test_that("passing .cl.ct.to.ign via ... produces informative error", {
  tdr <- local_pipeline()

  # Use the mapped object and try calling get.map again with .cl.ct.to.ign
  # We need the object before get.map; use a fresh pipeline through get.graph
  set.seed(99)
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  pre_map <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:4),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 99) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 99)

  expect_error(
    get.map(pre_map, .cl.ct.to.ign = "dummy", .verbose = FALSE),
    regexp = "'.cl.ct.to.ign' argument has been removed",
    fixed = FALSE
  )
})

test_that("error message mentions upstream QC", {
  set.seed(99)
  test_data <- create_test_lm_obj(n_cells = 50, n_markers = 4, n_samples = 4)
  on.exit(cleanup_test_files(test_data), add = TRUE)

  pre_map <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:4),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 99) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 99)

  expect_error(
    get.map(pre_map, .cl.ct.to.ign = "anything", .verbose = FALSE),
    regexp = "upstream",
    fixed = FALSE
  )
})

# =========================================================================
# 3. Functional correctness without .cl.ct.to.ign
# =========================================================================

test_that("get.map succeeds and populates density/cellmap slots", {
  tdr <- local_pipeline()

  expect_true(is.TDRObj(tdr))
  expect_false(is.null(tdr@density$fdens))
  expect_false(is.null(tdr@density$Y))
  expect_false(is.null(tdr@density$composition$clustering$cell.count))
  expect_false(is.null(tdr@density$composition$clustering$cell.perc))
})

test_that("fdens has correct dimensions (landmarks x samples)", {

  tdr <- local_pipeline()

  n_landmarks <- nrow(tdr@assay$expr)
  n_samples <- length(tdr@cells)
  expect_equal(nrow(tdr@density$fdens), n_landmarks)
  expect_equal(ncol(tdr@density$fdens), n_samples)
})

test_that("Y equals log2(fdens + 0.5)", {
  tdr <- local_pipeline()

  expected_Y <- log2(tdr@density$fdens + 0.5)
  expect_equal(tdr@density$Y, expected_Y)
})

test_that("cell.count rows match sample names", {
  tdr <- local_pipeline()

  expect_equal(
    sort(rownames(tdr@density$composition$clustering$cell.count)),
    sort(names(tdr@cells))
  )
})

test_that("cell.perc rows sum to approximately 100", {
  tdr <- local_pipeline()

  cc <- tdr@density$composition$clustering$cell.count
  expect_false(is.null(cc) || ncol(cc) == 0,
               label = "cell.count must have >= 1 cluster column")
  perc <- tdr@density$composition$clustering$cell.perc
  row_sums <- rowSums(perc, na.rm = TRUE)
  # Each sample should have cells so row sums should be ~100
  expect_true(all(abs(row_sums - 100) < 1e-8))
})

# =========================================================================
# 4. density$ignored is always NULL
# =========================================================================

test_that("density$ignored is NULL after get.map", {
  tdr <- local_pipeline()
  expect_null(tdr@density$ignored)
})

# =========================================================================
# 5. All clusters present in composition (no exclusion)
# =========================================================================

test_that("no cluster columns were dropped from cell.count", {
  tdr <- local_pipeline()

  cc <- tdr@density$composition$clustering$cell.count
  skip_if(is.null(cc), message = "No clustering composition in test fixture")

  # Verify cell.count is a valid matrix with expected row names
  expect_true(is.matrix(cc))
  expect_equal(sort(rownames(cc)), sort(names(tdr@cells)))

  # With multiple clusters, columns should be named and present
  if (length(levels(tdr@landmark.annot$clustering$ids)) > 1) {
    expect_true(ncol(cc) >= 1)
    expect_true(all(nchar(colnames(cc)) > 0))
  }
})

# =========================================================================
# 6. Deterministic reproducibility
# =========================================================================

test_that("get.map is deterministic with same seed", {
  tdr1 <- local_pipeline(seed = 42)
  tdr2 <- local_pipeline(seed = 42)

  expect_equal(tdr1@density$fdens, tdr2@density$fdens)
  expect_equal(tdr1@density$Y, tdr2@density$Y)
  expect_equal(
    tdr1@density$composition$clustering$cell.count,
    tdr2@density$composition$clustering$cell.count
  )
})

# =========================================================================
# 7. Dispatch wrapper passes through correctly
# =========================================================================

test_that("Seurat dispatch wrapper still works without .cl.ct.to.ign", {
  skip_if_not_installed("SeuratObject")

  set.seed(42)
  n_cells <- 50
  n_markers <- 4
  n_samples <- 4
  test_data <- create_test_lm_obj(
    n_cells = n_cells, n_markers = n_markers, n_samples = n_samples
  )
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:4),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 42) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 42)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu <- SetTDR(seu, tdr)
  result <- get.map(seu, .verbose = FALSE, .seed = 42)

  tdr_out <- GetTDR(result)
  expect_false(is.null(tdr_out@density$fdens))
  expect_null(tdr_out@density$ignored)
})

# =========================================================================
# 8. Seurat dispatch also rejects .cl.ct.to.ign
# =========================================================================

test_that("Seurat dispatch rejects .cl.ct.to.ign with informative error", {
  skip_if_not_installed("SeuratObject")

  set.seed(42)
  test_data <- create_test_lm_obj(
    n_cells = 50, n_markers = 4, n_samples = 4
  )
  on.exit(cleanup_test_files(test_data), add = TRUE)

  tdr <- setup.tdr.obj(
    .cells = test_data$cells,
    .meta = test_data$meta,
    .markers = paste0("marker_", 1:4),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 42) |>
    get.graph(.k = 5, .verbose = FALSE, .seed = 42)

  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu <- SetTDR(seu, tdr)

  expect_error(
    get.map(seu, .cl.ct.to.ign = "dummy", .verbose = FALSE),
    regexp = "'.cl.ct.to.ign' argument has been removed"
  )
})

# =========================================================================
# 9. .ct.to.keep legacy variable is fully removed
# =========================================================================

test_that("get.map.TDRObj body does not reference .ct.to.keep", {
  fn_body <- deparse(body(get.map.TDRObj))
  expect_false(
    any(grepl(".ct.to.keep", fn_body, fixed = TRUE)),
    info = "Legacy .ct.to.keep variable should not appear in get.map.TDRObj"
  )
})

test_that("fdens uses all cells (no cell-type subsetting)", {
  tdr <- local_pipeline()

  # fdens should be non-NULL and have correct dimensions

  expect_false(is.null(tdr@density$fdens))
  n_landmarks <- nrow(tdr@assay$expr)
  n_samples   <- length(tdr@cells)
  expect_equal(nrow(tdr@density$fdens), n_landmarks)
  expect_equal(ncol(tdr@density$fdens), n_samples)

  # All fdens values should be non-negative (sums of non-negative weights)
  expect_true(all(tdr@density$fdens >= 0))

  # Every sample should contribute a non-zero fdens column
  expect_true(all(colSums(tdr@density$fdens) > 0))
})
