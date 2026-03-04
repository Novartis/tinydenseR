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

# ── Shared fixture: build a cyto pipeline object and a full-cell label vector ──

.build_mode_b_fixture <- function(seed = 42, n_samples = 2, run_map = FALSE) {

  set.seed(seed)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  blob2 <- matrix(rnorm(n * m, mean = 5, sd = 0.5), nrow = n, ncol = m)

  sample_names <- paste0("sample", seq_len(n_samples))

  # Create expression matrices per sample and the full cell-label vector
  all_cell_ids    <- character(0)
  all_cell_labels <- character(0)

  .cells <- lapply(seq_len(n_samples), function(i) {
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    cell_ids <- paste0("s", i, "_", seq_len(nrow(mat)))
    dimnames(mat) <- list(cell_ids, paste0("M", seq_len(m)))

    # Assign labels: first half → "TypeA", second half → "TypeB"
    labels <- ifelse(seq_len(nrow(mat)) <= n, "TypeA", "TypeB")

    all_cell_ids    <<- c(all_cell_ids, cell_ids)
    all_cell_labels <<- c(all_cell_labels, labels)

    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    return(uri)
  })
  names(.cells) <- sample_names

  .meta <- data.frame(
    row.names = sample_names,
    group = rep(c("A", "B"), length.out = n_samples)
  )

  # Build the named character vector (Mode B input)
  cell_labels <- stats::setNames(all_cell_labels, all_cell_ids)

  obj <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("M", seq_len(m)),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(
      .k = 5,
      .scale = FALSE,
      .verbose = FALSE,
      .seed = seed
    )

  if (run_map) {
    obj <- get.map(obj, .verbose = FALSE, .seed = seed)
  }

  list(obj = obj, cell_labels = cell_labels)
}


# ──────────────────────────────────────────────────────────────────────
# U8: Mode dispatch: list → A, named char → B, other → error
# ──────────────────────────────────────────────────────────────────────

test_that("U8: mode dispatch selects correctly", {
  fix <- .build_mode_b_fixture()

  # unnamed numeric → error
  expect_error(
    celltyping(fix$obj, 1:10),
    "must be either"
  )

  # unnamed character → error
  expect_error(
    celltyping(fix$obj, c("A", "B")),
    "must be either"
  )

  # NULL → error
  expect_error(
    celltyping(fix$obj, NULL),
    "must be either"
  )
})

# ──────────────────────────────────────────────────────────────────────
# U2: Mode B: correct labels assigned to landmarks
# ──────────────────────────────────────────────────────────────────────

test_that("U2: Mode B assigns correct labels to landmarks", {
  fix <- .build_mode_b_fixture()
  obj <- fix$obj
  cell_labels <- fix$cell_labels

  obj2 <- celltyping(obj, cell_labels, .verbose = FALSE)

  ids <- obj2@landmark.annot$celltyping$ids
  expect_s3_class(ids, "factor")
  expect_true(all(ids %in% c("TypeA", "TypeB")))

  # Length must match number of landmarks
  expect_equal(length(ids), nrow(obj2@assay$expr))
})

# ──────────────────────────────────────────────────────────────────────
# U3: Mode B: non-landmark labels silently discarded
# ──────────────────────────────────────────────────────────────────────

test_that("U3: Mode B silently discards non-landmark labels", {
  fix <- .build_mode_b_fixture()
  obj <- fix$obj

  n_landmarks <- nrow(obj@assay$expr)
  n_total     <- sum(obj@config$sampling$n.cells)

  expect_true(n_landmarks < n_total)

  obj2 <- celltyping(obj, fix$cell_labels, .verbose = FALSE)
  expect_equal(length(obj2@landmark.annot$celltyping$ids), n_landmarks)
})

# ──────────────────────────────────────────────────────────────────────
# U4: Mode B: stop() when some landmarks have no matching name
# ──────────────────────────────────────────────────────────────────────

test_that("U4: Mode B errors when landmarks are not found in input", {
  fix <- .build_mode_b_fixture()
  obj <- fix$obj

  # Remove some cells from the labels vector
  bad_labels <- fix$cell_labels[1:(length(fix$cell_labels) - 5)]
  # Need correct length but wrong names
  # Actually just remove names of cells that are landmarks
  lm_names <- rownames(obj@assay$expr)
  sample_prefixes <- paste0("^", unique(names(obj@config$key)), "_")
  original_ids <- sub(
    pattern = paste0(sample_prefixes, collapse = "|"),
    replacement = "",
    x = lm_names
  )
  # Remove some landmark original IDs from the labels
  drop_ids <- head(original_ids, 3)
  bad_labels <- fix$cell_labels
  # Rename cells so landmarks can't be found
  names(bad_labels)[names(bad_labels) %in% drop_ids] <-
    paste0("FAKE_", drop_ids)

  expect_error(
    celltyping(obj, bad_labels, .verbose = FALSE),
    "Could not find labels"
  )
})

# ──────────────────────────────────────────────────────────────────────
# U5: Mode B: stop() when length doesn't match total cell count
# ──────────────────────────────────────────────────────────────────────

test_that("U5: Mode B errors on length mismatch", {
  fix <- .build_mode_b_fixture()

  too_short <- fix$cell_labels[1:10]
  expect_error(
    celltyping(fix$obj, too_short, .verbose = FALSE),
    "does not match total cell count"
  )
})

# ──────────────────────────────────────────────────────────────────────
# U6: Mode B: stop() when names have duplicates
# ──────────────────────────────────────────────────────────────────────

test_that("U6: Mode B errors on duplicate cell IDs", {
  fix <- .build_mode_b_fixture()

  dup_labels <- fix$cell_labels
  names(dup_labels)[2] <- names(dup_labels)[1]  # introduce duplicate

  expect_error(
    celltyping(fix$obj, dup_labels, .verbose = FALSE),
    "Duplicate cell IDs"
  )
})

# ──────────────────────────────────────────────────────────────────────
# U7: Mode B: stop() when vector is unnamed
# ──────────────────────────────────────────────────────────────────────

test_that("U7: unnamed character vector gives clear error", {
  fix <- .build_mode_b_fixture()

  unnamed <- unname(fix$cell_labels)
  expect_error(
    celltyping(fix$obj, unnamed, .verbose = FALSE),
    "must be either"
  )
})

# ──────────────────────────────────────────────────────────────────────
# U9: $mode slot is correct for both modes
# ──────────────────────────────────────────────────────────────────────

test_that("U9: $mode slot is cluster_map for A, cell_labels for B", {
  fix <- .build_mode_b_fixture()
  obj <- fix$obj

  # Mode A
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  half <- max(1, floor(n / 2))
  ct_map <- list(
    "TypeA" = cls[seq_len(half)],
    "TypeB" = cls[(half + 1):n]
  )
  obj_a <- celltyping(obj, ct_map, .verbose = FALSE)
  expect_identical(obj_a@landmark.annot$celltyping$mode, "cluster_map")
  expect_true(!is.null(obj_a@landmark.annot$celltyping$map))

  # Mode B
  obj_b <- celltyping(obj, fix$cell_labels, .verbose = FALSE)
  expect_identical(obj_b@landmark.annot$celltyping$mode, "cell_labels")
  expect_null(obj_b@landmark.annot$celltyping$map)
})

# ──────────────────────────────────────────────────────────────────────
# I1: Full pipeline Mode B with get.map() — verify downstream slots
# ──────────────────────────────────────────────────────────────────────

test_that("I1: Mode B after get.map refreshes all downstream slots", {
  fix <- .build_mode_b_fixture(run_map = TRUE)
  obj <- fix$obj

  obj2 <- celltyping(obj, fix$cell_labels, .verbose = FALSE)

  # Cell-level IDs should be populated
  ct_ids <- .tdr_get_map_slot_all(obj2, "celltyping.ids")
  expect_true(!is.null(ct_ids))
  all_labels <- unique(unlist(lapply(ct_ids, unique)))
  expect_true(any(c("TypeA", "TypeB") %in% all_labels))

  # Composition matrices
  cell_count <- obj2@density$composition$celltyping$cell.count
  cell_perc  <- obj2@density$composition$celltyping$cell.perc
  expect_true(is.matrix(cell_count))
  expect_true(is.matrix(cell_perc))
  expect_true(any(c("TypeA", "TypeB") %in% colnames(cell_count)))

  # Median expression + pheatmap
  expect_true(is.matrix(obj2@results$celltyping$median.exprs))
  expect_s3_class(obj2@results$celltyping$pheatmap, "pheatmap")
})

# ──────────────────────────────────────────────────────────────────────
# I4: Switch from Mode A to Mode B on same object
# ──────────────────────────────────────────────────────────────────────

test_that("I4: switching Mode A → Mode B fully replaces slots", {
  fix <- .build_mode_b_fixture(run_map = TRUE)
  obj <- fix$obj

  # First apply Mode A
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  half <- max(1, floor(n / 2))
  ct_map <- list("Alpha" = cls[seq_len(half)], "Beta" = cls[(half + 1):n])
  obj_a <- celltyping(obj, ct_map, .verbose = FALSE)
  expect_identical(obj_a@landmark.annot$celltyping$mode, "cluster_map")
  expect_true(all(obj_a@landmark.annot$celltyping$ids %in% c("Alpha", "Beta")))

  # Switch to Mode B
  obj_b <- celltyping(obj_a, fix$cell_labels, .verbose = FALSE)
  expect_identical(obj_b@landmark.annot$celltyping$mode, "cell_labels")
  expect_null(obj_b@landmark.annot$celltyping$map)
  expect_true(all(obj_b@landmark.annot$celltyping$ids %in% c("TypeA", "TypeB")))

  # Median expression rownames should reflect Mode B labels
  expect_true(all(
    sort(rownames(obj_b@results$celltyping$median.exprs)) ==
    sort(c("TypeA", "TypeB"))
  ))
})

# ──────────────────────────────────────────────────────────────────────
# Mode B warning when all landmarks get same label
# ──────────────────────────────────────────────────────────────────────

test_that("Mode B warns when all landmarks receive the same label", {
  fix <- .build_mode_b_fixture()
  obj <- fix$obj

  # Make all labels the same
  uniform_labels <- fix$cell_labels
  uniform_labels[] <- "OnlyType"

  expect_warning(
    celltyping(obj, uniform_labels, .verbose = FALSE),
    "same cell-type label"
  )
})
