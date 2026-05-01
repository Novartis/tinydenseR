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
# Tests: get.map() with .ref.obj preserves existing named celltyping solutions
#
# The Symphony code path in get.map() now preserves pre-existing named
# solutions in @landmark.annot$celltyping instead of wiping the entire
# list.  These tests validate that invariant.
# =========================================================================

library(testthat)
library(tinydenseR)

# ── Shared fixture: minimal cyto pipeline through get.map() + celltyping ──
.build_typed_obj <- function(seed = 42, n_samples = 4) {
  set.seed(seed)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  blob2 <- matrix(rnorm(n * m, mean = 5, sd = 0.5), nrow = n, ncol = m)

  .cells <- lapply(seq_len(n_samples), function(i) {
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    dimnames(mat) <- list(
      paste0("s", i, "_", seq_len(nrow(mat))),
      paste0("M", seq_len(m))
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    return(uri)
  })
  names(.cells) <- paste0("sample", seq_len(n_samples))

  .meta <- data.frame(
    row.names = names(.cells),
    group = rep(c("A", "B"), length.out = n_samples)
  )

  obj <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("M", seq_len(m)),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = seed) |>
    get.map(.verbose = FALSE, .seed = seed, .cache.on.disk = FALSE)

  # Add a manual celltyping solution
  cls <- levels(obj@landmark.annot$clustering$ids)
  half <- floor(length(cls) / 2)
  ct_map <- list(
    "TypeA" = cls[seq_len(half)],
    "TypeB" = cls[(half + 1):length(cls)]
  )
  obj <- celltyping(obj, ct_map, .verbose = FALSE, .name = "manual_v1")
  obj
}

# =========================================================================
# 1. Pre-loop preservation: clearing $ids does NOT destroy named solutions
# =========================================================================

test_that("clearing $ids preserves named celltyping solutions", {
  obj <- .build_typed_obj()

  # Verify the named solution exists before the simulated clear

  solutions_before <- list_celltyping_solutions(obj)
  expect_true("manual_v1" %in% solutions_before)
  expect_false(is.null(obj@landmark.annot$celltyping$ids))

  # Simulate the pre-loop logic from get.map() when .ref.obj is provided:
  #   if (is.null(.tdr.obj@landmark.annot$celltyping)) {
  #     .tdr.obj@landmark.annot$celltyping <- list()
  #   }
  #   .tdr.obj@landmark.annot$celltyping$ids <- NULL
  if (is.null(obj@landmark.annot$celltyping)) {
    obj@landmark.annot$celltyping <- list()
  }
  obj@landmark.annot$celltyping$ids <- NULL

  # Named solutions survive; $ids is gone
  expect_null(obj@landmark.annot$celltyping$ids)
  expect_true("manual_v1" %in% names(obj@landmark.annot$celltyping))
  expect_true(is.factor(obj@landmark.annot$celltyping[["manual_v1"]]))
  expect_equal(
    length(obj@landmark.annot$celltyping[["manual_v1"]]),
    nrow(obj@assay$expr)
  )
})

# =========================================================================
# 2. Post-loop rebuild: new $ids + named solution coexist with old solutions
# =========================================================================

test_that("post-loop rebuild adds Symphony solution alongside existing ones", {
  obj <- .build_typed_obj()

  saved_manual <- obj@landmark.annot$celltyping[["manual_v1"]]

  # Simulate the pre-loop clear
  obj@landmark.annot$celltyping$ids <- NULL

  # Simulate the post-loop rebuild (as if Symphony assigned labels)
  n_lm <- nrow(obj@assay$expr)
  fake_symphony_labels <- factor(
    sample(c("Mono", "T_cell", "B_cell"), n_lm, replace = TRUE),
    levels = c("B_cell", "Mono", "T_cell")
  )
  names(fake_symphony_labels) <- rownames(obj@assay$expr)

  obj@landmark.annot$celltyping$ids <- fake_symphony_labels
  obj@landmark.annot$celltyping[["cell_type_L1"]] <- fake_symphony_labels

  # Both solutions coexist
  solutions <- list_celltyping_solutions(obj)
  expect_true("manual_v1" %in% solutions)
  expect_true("cell_type_L1" %in% solutions)

  # Active solution is the Symphony one
  expect_identical(obj@landmark.annot$celltyping$ids, fake_symphony_labels)

  # Old solution is untouched
  expect_identical(obj@landmark.annot$celltyping[["manual_v1"]], saved_manual)
})

# =========================================================================
# 3. set_active_celltyping() can restore original solution after ref override
# =========================================================================

test_that("set_active_celltyping restores manual solution after simulated ref override", {
  obj <- .build_typed_obj()

  saved_manual <- obj@landmark.annot$celltyping[["manual_v1"]]

  # Simulate ref-based override
  obj@landmark.annot$celltyping$ids <- NULL
  n_lm <- nrow(obj@assay$expr)
  fake_labels <- factor(
    sample(c("DC", "NK"), n_lm, replace = TRUE),
    levels = c("DC", "NK")
  )
  names(fake_labels) <- rownames(obj@assay$expr)
  obj@landmark.annot$celltyping$ids <- fake_labels
  obj@landmark.annot$celltyping[["ref_ct"]] <- fake_labels

  # Now switch back to manual
  obj2 <- set_active_celltyping(obj, "manual_v1", .verbose = FALSE)

  expect_identical(obj2@landmark.annot$celltyping$ids, saved_manual)
  # Both solutions still listed
  expect_true(all(c("manual_v1", "ref_ct") %in% list_celltyping_solutions(obj2)))
})

# =========================================================================
# 4. Preservation from NULL: when no prior celltyping exists
# =========================================================================

test_that("preservation logic handles NULL celltyping gracefully", {
  obj <- .build_typed_obj()

  # Wipe celltyping entirely (simulates object that never had celltyping)
  obj@landmark.annot$celltyping <- NULL

  # Simulate the pre-loop logic
  if (is.null(obj@landmark.annot$celltyping)) {
    obj@landmark.annot$celltyping <- list()
  }
  obj@landmark.annot$celltyping$ids <- NULL

  # Should be an empty list with no $ids
  expect_type(obj@landmark.annot$celltyping, "list")
  expect_equal(length(obj@landmark.annot$celltyping), 0L)
})

# =========================================================================
# 5. Warning message updated (no longer says "Removing")
# =========================================================================

test_that("get.map .ref.obj warning mentions preserving solutions", {
  obj <- .build_typed_obj()

  # Provide a syntactically valid but incomplete ref.obj.
  # This will pass the initial checks (is list, has Z_corr, has column)
  # but will error downstream (e.g. during annoy_build or mapQuery).
  # We catch the warning before the error.
  fake_ref <- list(
    Z_corr = matrix(rnorm(20), nrow = 2, ncol = 10),
    meta_data = data.frame(cell_type = rep("A", 10))
  )

  # Cytometry objects reject .ref.obj, so we temporarily fake RNA assay type
  # to reach the warning. The downstream error is expected.
  obj@config$assay.type <- "RNA"

  expect_warning(
    tryCatch(
      get.map(obj, .ref.obj = fake_ref, .verbose = FALSE),
      error = function(e) NULL
    ),
    "Previous named solutions are preserved"
  )
})
