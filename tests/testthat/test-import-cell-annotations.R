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

# ======================================================================
# test-import-cell-annotations.R
#
# Comprehensive test suite for automatic multi-column cell-level
# annotation ingestion via import_cell_annotations().
# Covers:
#   - Column selection criteria (char/factor, not numeric, not
#     .sample.var, not "ids", not ID-like, not single-value,
#     cell-level not sample-level)
#   - Landmark storage (each column → named solution in
#     @landmark.annot$celltyping$<colname>)
#   - Active annotation (.celltype.vec → active $ids; NULL → $ids
#     unchanged)
#   - Full pipeline integration (import → get.map → cell-level IDs)
#   - Switching between imported solutions
#   - Backward compatibility (manual celltyping Mode A / Mode B)
#   - list_celltyping_solutions() utility
#   - Edge cases (no categorical columns, all-NA, single-sample,
#     pre-map object, column name collision)
# ======================================================================

library(testthat)
library(tinydenseR)

# ======================================================================
# FIXTURE BUILDERS
# ======================================================================

# Build a graph-stage object with full cell metadata available
# Returns list(obj = TDRObj, cell_meta = data.frame)
.build_graph_with_meta <- function(seed = 42, n_samples = 4,
                                   n_per_sample = 200,
                                   n_markers = 5) {
  set.seed(seed)
  blob1 <- matrix(rnorm(n_per_sample * n_markers, mean = 0, sd = 0.5),
                  nrow = n_per_sample, ncol = n_markers)
  blob2 <- matrix(rnorm(n_per_sample * n_markers, mean = 5, sd = 0.5),
                  nrow = n_per_sample, ncol = n_markers)

  sample_names <- paste0("sample", seq_len(n_samples))
  all_cell_ids <- character(0)
  cell_sample  <- character(0)

  .cells <- lapply(seq_len(n_samples), function(i) {
    sn <- sample_names[i]
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    cids <- paste0("s", i, "_cell_", seq_len(nrow(mat)))
    dimnames(mat) <- list(
      cids,
      paste0("M", seq_len(n_markers))
    )
    all_cell_ids <<- c(all_cell_ids, cids)
    cell_sample  <<- c(cell_sample, rep(sn, nrow(mat)))
    uri <- tempfile(fileext = ".RDS")
    saveRDS(mat, file = uri, compress = FALSE)
    uri
  })
  names(.cells) <- sample_names

  .meta <- data.frame(
    row.names = sample_names,
    group = rep(c("A", "B"), length.out = n_samples)
  )

  obj <- setup.tdr.obj(
    .cells = .cells, .meta = .meta,
    .markers = paste0("M", seq_len(n_markers)),
    .assay.type = "cyto", .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = seed)

  n_total <- n_samples * n_per_sample * 2  # blob1 + blob2

  # Build cell-level metadata with various column types
  cell_meta <- data.frame(
    row.names       = all_cell_ids,
    sample_id       = cell_sample,
    # Cell-level character: 2 labels, varies within samples
    cell_type_l1    = rep(c("TypeA", "TypeB"), length.out = n_total),
    # Cell-level factor: 3 labels, varies within samples
    cell_type_l2    = factor(rep(c("Sub1", "Sub2", "Sub3"),
                                 length.out = n_total)),
    # Numeric column (should NOT be imported)
    nCount_RNA      = runif(n_total, 1000, 5000),
    # Logical column (should NOT be imported)
    is_doublet      = sample(c(TRUE, FALSE), n_total, replace = TRUE),
    # Integer column (should NOT be imported)
    n_genes         = sample(500:2000, n_total, replace = TRUE),
    # Sample-level constant column (should NOT be imported)
    batch           = rep(c("batch1", "batch2"),
                          each = n_per_sample * 2,
                          length.out = n_total),
    # ID-like column: all unique (should NOT be imported)
    barcode         = paste0("BC_", seq_len(n_total)),
    # Single-value column (should NOT be imported)
    species         = rep("human", n_total),
    stringsAsFactors = FALSE
  )
  # batch is constant within each sample (sample-level), confirm:
  # sample1 → batch1, sample2 → batch2, etc.
  cell_meta$batch <- rep(
    rep(c("batch1", "batch2"), length.out = n_samples),
    each = n_per_sample * 2
  )

  list(obj = obj, cell_meta = cell_meta)
}

# Build a cluster→celltype map from whatever clusters exist
.make_ct_map <- function(obj) {
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  stopifnot(n >= 2)
  half <- floor(n / 2)
  list(
    "TypeA" = cls[seq_len(half)],
    "TypeB" = cls[(half + 1):n]
  )
}


# ======================================================================
# CATEGORY 1 — COLUMN SELECTION CRITERIA
# ======================================================================

test_that("ICA1.1: character and factor columns are imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_true("cell_type_l1" %in% solutions)
  expect_true("cell_type_l2" %in% solutions)
})

test_that("ICA1.2: numeric columns are NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("nCount_RNA" %in% solutions)
  expect_false("n_genes" %in% solutions)
})

test_that("ICA1.3: logical columns are NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("is_doublet" %in% solutions)
})

test_that("ICA1.4: .sample.var column is NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("sample_id" %in% solutions)
})

test_that("ICA1.5: ID-like columns (all unique values) are NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("barcode" %in% solutions)
})

test_that("ICA1.6: single-value columns are NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("species" %in% solutions)
})

test_that("ICA1.7: column named 'ids' is excluded", {
  fix <- .build_graph_with_meta()
  fix$cell_meta$ids <- rep(c("a", "b"), length.out = nrow(fix$cell_meta))
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("ids" %in% solutions)
})

test_that("ICA1.8: sample-level constant columns are NOT imported", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  # batch is constant within each sample → sample-level, not imported
  expect_false("batch" %in% solutions)
})


# ======================================================================
# CATEGORY 2 — LANDMARK STORAGE
# ======================================================================

test_that("ICA2.1: each imported column stored as named solution", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  # Each qualifying column is a named entry in @landmark.annot$celltyping
  expect_true(!is.null(obj@landmark.annot$celltyping$cell_type_l1))
  expect_true(!is.null(obj@landmark.annot$celltyping$cell_type_l2))
})

test_that("ICA2.2: stored solutions are factors of correct length", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)

  n_lm <- nrow(obj@assay$expr)
  sol_l1 <- obj@landmark.annot$celltyping$cell_type_l1
  sol_l2 <- obj@landmark.annot$celltyping$cell_type_l2

  expect_true(is.factor(sol_l1))
  expect_true(is.factor(sol_l2))
  expect_equal(length(sol_l1), n_lm)
  expect_equal(length(sol_l2), n_lm)
})

test_that("ICA2.3: landmark labels come from Mode B (original cell's label)", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)

  sol <- obj@landmark.annot$celltyping$cell_type_l1
  # All labels must be one of the original values
  expect_true(all(as.character(sol) %in% c("TypeA", "TypeB")))
})


# ======================================================================
# CATEGORY 3 — ACTIVE ANNOTATION
# ======================================================================

test_that("ICA3.1: .celltype.vec activates the specified column as $ids", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)

  expect_identical(
    obj@landmark.annot$celltyping$ids,
    obj@landmark.annot$celltyping$cell_type_l1
  )
})

test_that("ICA3.2: without .celltype.vec, $ids is unchanged", {
  fix <- .build_graph_with_meta()
  ids_before <- fix$obj@landmark.annot$celltyping$ids

  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = NULL,
                                 .verbose = FALSE)

  expect_identical(obj@landmark.annot$celltyping$ids, ids_before)
})

test_that("ICA3.3: .celltype.vec not among imported columns → $ids unchanged", {
  fix <- .build_graph_with_meta()
  ids_before <- fix$obj@landmark.annot$celltyping$ids

  # "nCount_RNA" is numeric → won't be imported
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "nCount_RNA",
                                 .verbose = FALSE)

  expect_identical(obj@landmark.annot$celltyping$ids, ids_before)
})


# ======================================================================
# CATEGORY 4 — FULL PIPELINE INTEGRATION
# ======================================================================

test_that("ICA4.1: import → get.map → cell-level IDs from active annotation", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  obj <- get.map(obj, .verbose = FALSE, .seed = 42)

  # Cell-level celltyping IDs should exist
  for (sn in names(obj@cells)) {
    ct_ids <- tinydenseR:::.tdr_get_map_slot(obj, "celltyping", sn)
    expect_true(!is.null(ct_ids))
    # All labels from the active vocabulary
    expect_true(all(ct_ids %in% c(levels(obj@landmark.annot$celltyping$ids), NA)))
  }
})

test_that("ICA4.2: switching imported solutions updates cell-level IDs", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  obj <- get.map(obj, .verbose = FALSE, .seed = 42)

  # Switch to cell_type_l2

  obj2 <- set_active_celltyping(obj, "cell_type_l2", .verbose = FALSE)

  # Active ids should now reflect l2
  expect_identical(
    obj2@landmark.annot$celltyping$ids,
    obj2@landmark.annot$celltyping$cell_type_l2
  )

  # Composition should reflect new labels
  comp <- obj2@density$composition$celltyping$cell.count
  expect_true(is.matrix(comp))
  expect_true(all(
    levels(obj2@landmark.annot$celltyping$ids) %in% colnames(comp)
  ))
})


# ======================================================================
# CATEGORY 5 — SWITCHING BETWEEN SOLUTIONS
# ======================================================================

test_that("ICA5.1: switch between imported solutions is lossless", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  obj <- get.map(obj, .verbose = FALSE, .seed = 42)

  sol_l1 <- obj@landmark.annot$celltyping$cell_type_l1
  sol_l2 <- obj@landmark.annot$celltyping$cell_type_l2

  # Switch l1 → l2 → l1
  obj2 <- set_active_celltyping(obj, "cell_type_l2", .verbose = FALSE)
  obj3 <- set_active_celltyping(obj2, "cell_type_l1", .verbose = FALSE)

  expect_identical(obj3@landmark.annot$celltyping$cell_type_l1, sol_l1)
  expect_identical(obj3@landmark.annot$celltyping$cell_type_l2, sol_l2)
  expect_identical(obj3@landmark.annot$celltyping$ids, sol_l1)
})

test_that("ICA5.2: composition matrices update on switch", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  obj <- get.map(obj, .verbose = FALSE, .seed = 42)

  comp_l1 <- obj@density$composition$celltyping$cell.count
  obj2 <- set_active_celltyping(obj, "cell_type_l2", .verbose = FALSE)
  comp_l2 <- obj2@density$composition$celltyping$cell.count

  # Different label sets → different column names
  expect_true(all(c("TypeA", "TypeB") %in% colnames(comp_l1)))
  expect_true(all(c("Sub1", "Sub2", "Sub3") %in% colnames(comp_l2)))
})


# ======================================================================
# CATEGORY 6 — BACKWARD COMPATIBILITY
# ======================================================================

test_that("ICA6.1: manual celltyping Mode A still works after import", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "manual_v1", .verbose = FALSE)

  expect_true("manual_v1" %in% list_celltyping_solutions(obj2))
  # Prior imported solutions still present
  expect_true("cell_type_l1" %in% list_celltyping_solutions(obj2))
  expect_true("cell_type_l2" %in% list_celltyping_solutions(obj2))
})

test_that("ICA6.2: manual celltyping Mode B still works after import", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)

  # Build a Mode B named character vector from cell_meta
  vec <- stats::setNames(
    as.character(fix$cell_meta$cell_type_l1),
    rownames(fix$cell_meta)
  )
  obj2 <- celltyping(obj, .celltyping.map = vec,
                     .name = "manual_b", .verbose = FALSE)

  expect_true("manual_b" %in% list_celltyping_solutions(obj2))
  expect_true(is.factor(obj2@landmark.annot$celltyping$ids))
})

test_that("ICA6.3: import_cell_annotations is a no-op with no categorical cols", {
  fix <- .build_graph_with_meta()
  # Build metadata with only numeric and sample-level columns
  meta_numeric <- data.frame(
    row.names = rownames(fix$cell_meta),
    sample_id = fix$cell_meta$sample_id,
    score1 = runif(nrow(fix$cell_meta)),
    score2 = rnorm(nrow(fix$cell_meta)),
    stringsAsFactors = FALSE
  )

  ids_before <- fix$obj@landmark.annot$celltyping$ids
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = meta_numeric,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)

  # No solutions should have been added
  expect_equal(length(list_celltyping_solutions(obj)), 0L)
  expect_identical(obj@landmark.annot$celltyping$ids, ids_before)
})


# ======================================================================
# CATEGORY 7 — EDGE CASES
# ======================================================================

test_that("ICA7.1: no qualifying columns returns object unchanged", {
  fix <- .build_graph_with_meta()
  # Only sample_id (excluded) and numeric columns
  meta_only_numeric <- data.frame(
    row.names = rownames(fix$cell_meta),
    sample_id = fix$cell_meta$sample_id,
    x = runif(nrow(fix$cell_meta)),
    stringsAsFactors = FALSE
  )

  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = meta_only_numeric,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  expect_equal(length(list_celltyping_solutions(obj)), 0L)
})

test_that("ICA7.2: column with all NA values (single non-NA unique) skipped", {
  fix <- .build_graph_with_meta()
  fix$cell_meta$all_na <- NA_character_
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("all_na" %in% solutions)
})

test_that("ICA7.3: single-sample object works", {
  fix <- .build_graph_with_meta(n_samples = 1)
  # With 1 sample, cell_type_l1 varies within that sample
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  expect_true("cell_type_l1" %in% list_celltyping_solutions(obj))
  expect_identical(
    obj@landmark.annot$celltyping$ids,
    obj@landmark.annot$celltyping$cell_type_l1
  )
})

test_that("ICA7.4: works on object where get.map() hasn't run", {
  fix <- .build_graph_with_meta()
  # obj is at graph stage, no get.map yet
  expect_null(fix$obj@density$norm)

  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)

  # Solutions stored correctly
  expect_true("cell_type_l1" %in% list_celltyping_solutions(obj))
  # No crash — .refresh_celltyping is a no-op since norm is NULL
  expect_null(obj@density$norm)
})

test_that("ICA7.5: column name collision with existing named solution", {
  fix <- .build_graph_with_meta()
  # Manually add a solution named "cell_type_l1"
  ct_map <- .make_ct_map(fix$obj)
  fix$obj <- celltyping(fix$obj, ct_map, .name = "cell_type_l1",
                        .verbose = FALSE)

  # import_cell_annotations should overwrite the named solution
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)

  # The solution should reflect Mode B (cell_type_l1 labels), not Mode A
  sol <- obj@landmark.annot$celltyping$cell_type_l1
  expect_true(all(as.character(sol) %in% c("TypeA", "TypeB")))
})

test_that("ICA7.6: verbose messages about imported columns", {
  fix <- .build_graph_with_meta()
  expect_message(
    import_cell_annotations(fix$obj,
                            .cell.meta = fix$cell_meta,
                            .sample.var = "sample_id",
                            .verbose = TRUE),
    "importing.*column"
  )
})

test_that("ICA7.7: verbose message when no qualifying columns", {
  fix <- .build_graph_with_meta()
  meta_empty <- data.frame(
    row.names = rownames(fix$cell_meta),
    sample_id = fix$cell_meta$sample_id,
    stringsAsFactors = FALSE
  )
  expect_message(
    import_cell_annotations(fix$obj,
                            .cell.meta = meta_empty,
                            .sample.var = "sample_id",
                            .verbose = TRUE),
    "no.*categorical"
  )
})


# ======================================================================
# CATEGORY 8 — list_celltyping_solutions()
# ======================================================================

test_that("ICA8.1: returns correct names after import", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_true("cell_type_l1" %in% solutions)
  expect_true("cell_type_l2" %in% solutions)
})

test_that("ICA8.2: returns empty character(0) when no solutions stored", {
  fix <- .build_graph_with_meta()
  solutions <- list_celltyping_solutions(fix$obj)
  expect_equal(length(solutions), 0L)
  expect_identical(solutions, character(0))
})

test_that("ICA8.3: does not include 'ids' in returned names", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .celltype.vec = "cell_type_l1",
                                 .verbose = FALSE)
  solutions <- list_celltyping_solutions(obj)
  expect_false("ids" %in% solutions)
})

test_that("ICA8.4: includes both imported and manually added solutions", {
  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .name = "manual_v1", .verbose = FALSE)

  solutions <- list_celltyping_solutions(obj2)
  expect_true("cell_type_l1" %in% solutions)
  expect_true("cell_type_l2" %in% solutions)
  expect_true("manual_v1" %in% solutions)
})


# ======================================================================
# CATEGORY 9 — INPUT VALIDATION
# ======================================================================

test_that("ICA9.1: errors when .cell.meta is not a data.frame", {
  fix <- .build_graph_with_meta()
  expect_error(
    import_cell_annotations(fix$obj,
                            .cell.meta = "not_a_df",
                            .sample.var = "sample_id",
                            .verbose = FALSE),
    "data\\.frame"
  )
})

test_that("ICA9.2: errors when .sample.var not in .cell.meta", {
  fix <- .build_graph_with_meta()
  expect_error(
    import_cell_annotations(fix$obj,
                            .cell.meta = fix$cell_meta,
                            .sample.var = "nonexistent",
                            .verbose = FALSE),
    "not found"
  )
})

test_that("ICA9.3: errors when .cell.meta has no rownames", {
  fix <- .build_graph_with_meta()
  meta_no_rn <- fix$cell_meta
  rownames(meta_no_rn) <- NULL
  expect_error(
    import_cell_annotations(fix$obj,
                            .cell.meta = meta_no_rn,
                            .sample.var = "sample_id",
                            .verbose = FALSE),
    "rownames"
  )
})


# ======================================================================
# CATEGORY 10 — DISPATCH WRAPPERS
# ======================================================================

test_that("ICA10.1: import_cell_annotations dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")

  fix <- .build_graph_with_meta()
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat <- SetTDR(seurat, fix$obj)

  result <- suppressWarnings(
    import_cell_annotations(seurat,
                            .cell.meta = fix$cell_meta,
                            .sample.var = "sample_id",
                            .celltype.vec = "cell_type_l1",
                            .verbose = FALSE)
  )
  tdr_out <- GetTDR(result)
  expect_true("cell_type_l1" %in% list_celltyping_solutions(tdr_out))
})

test_that("ICA10.2: list_celltyping_solutions dispatches on Seurat", {
  skip_if_not_installed("SeuratObject")

  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 1, 1),
    dims = c(3, 3),
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))
  )
  seurat <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat <- SetTDR(seurat, obj)

  solutions <- list_celltyping_solutions(seurat)
  expect_true("cell_type_l1" %in% solutions)
})

test_that("ICA10.3: import_cell_annotations dispatches on SCE", {
  skip_if_not_installed("SingleCellExperiment")

  fix <- .build_graph_with_meta()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 3))
  )
  sce <- SetTDR(sce, fix$obj)

  result <- import_cell_annotations(sce,
                                    .cell.meta = fix$cell_meta,
                                    .sample.var = "sample_id",
                                    .celltype.vec = "cell_type_l2",
                                    .verbose = FALSE)
  tdr_out <- GetTDR(result)
  expect_true("cell_type_l2" %in% list_celltyping_solutions(tdr_out))
})

test_that("ICA10.4: list_celltyping_solutions dispatches on SCE", {
  skip_if_not_installed("SingleCellExperiment")

  fix <- .build_graph_with_meta()
  obj <- import_cell_annotations(fix$obj,
                                 .cell.meta = fix$cell_meta,
                                 .sample.var = "sample_id",
                                 .verbose = FALSE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(1, nrow = 3, ncol = 3))
  )
  sce <- SetTDR(sce, obj)

  solutions <- list_celltyping_solutions(sce)
  expect_true("cell_type_l1" %in% solutions)
})
