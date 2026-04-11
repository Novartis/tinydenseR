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

# ──────────────────────────────────────────────────────────────────────
# Tests for unified cellmap structure (schema v3)
#
# Validates:
# - Path string attribute preservation
# - Transparent lazy-load
# - Round-trip all 4 slots
# - Backward compat: old slot names
# - Migration: old @density$.cache → @cellmap path strings
# - Validate/info/cleanup with new structure
# ──────────────────────────────────────────────────────────────────────

# ── Path string attribute preservation ──────────────────────────────

test_that("path string attributes survive storage in @cellmap list", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_attr")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write(c(a = "cl.01"), cache_dir, "clustering", "s1")

  # Store in a cellmap-like list structure
  cm <- list(clustering = list(ids = list(s1 = p1)))

  # Verify attribute survives

  expect_identical(attr(cm$clustering$ids$s1, "schema_v"), .TDR_CACHE_SCHEMA_VERSION)
  expect_true(attr(cm$clustering$ids$s1, "bytes") > 0)
})

test_that("path string attributes survive TDRObj construction", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_attr_obj")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write(c(a = "cl.01"), cache_dir, "clustering", "s1")

  obj <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  val <- obj@cellmap$clustering$ids$s1
  expect_identical(attr(val, "schema_v"), .TDR_CACHE_SCHEMA_VERSION)
  expect_true(attr(val, "bytes") > 0)
})

# ── Transparent lazy-load ────────────────────────────────────────────

test_that("lazy-load: .tdr_get_map_slot transparently reads on-disk data", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_lazy")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  data_cl <- c(a = "cluster.01", b = "cluster.02")
  data_fg <- Matrix::rsparsematrix(10, 5, 0.3)
  data_nl <- matrix(1:6, nrow = 3, ncol = 2)
  data_ct <- c(a = "TypeA", b = "TypeB")

  p_cl <- .tdr_cache_write(data_cl, cache_dir, "clustering",   "s1")
  p_fg <- .tdr_cache_write(data_fg, cache_dir, "fuzzy.graphs", "s1")
  p_nl <- .tdr_cache_write(data_nl, cache_dir, "nearest.lm",   "s1")
  p_ct <- .tdr_cache_write(data_ct, cache_dir, "celltyping",   "s1")

  obj <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      celltyping   = list(ids = list(s1 = p_ct)),
      nearest.lm   = list(s1 = p_nl),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  expect_identical(.tdr_get_map_slot(obj, "clustering",   "s1"), data_cl)
  expect_identical(.tdr_get_map_slot(obj, "celltyping",   "s1"), data_ct)
  expect_identical(.tdr_get_map_slot(obj, "nearest.lm",   "s1"), data_nl)
  expect_true(Matrix::all.equal(
    .tdr_get_map_slot(obj, "fuzzy.graphs", "s1"), data_fg
  ))
})

# ── Round-trip all 4 slots ───────────────────────────────────────────

test_that("round-trip: write all 4 slot types, read back identically", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_roundtrip")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  samples <- c("s1", "s2")
  data <- list(
    clustering   = list(s1 = c(a = "1", b = "2"), s2 = c(c = "1")),
    celltyping   = list(s1 = c(a = "A", b = "B"), s2 = c(c = "A")),
    nearest.lm   = list(s1 = matrix(1:4, 2, 2), s2 = matrix(5:6, 1, 2)),
    fuzzy.graphs = list(
      s1 = Matrix::rsparsematrix(10, 5, 0.3),
      s2 = Matrix::rsparsematrix(8, 5, 0.2)
    )
  )

  cm <- list(
    clustering = list(ids = list()),
    celltyping = list(ids = list()),
    nearest.lm   = list(),
    fuzzy.graphs = list()
  )

  for (sn in samples) {
    cm$clustering$ids[[sn]] <- .tdr_cache_write(data$clustering[[sn]], cache_dir, "clustering", sn)
    cm$celltyping$ids[[sn]] <- .tdr_cache_write(data$celltyping[[sn]], cache_dir, "celltyping", sn)
    cm$nearest.lm[[sn]]     <- .tdr_cache_write(data$nearest.lm[[sn]], cache_dir, "nearest.lm", sn)
    cm$fuzzy.graphs[[sn]]   <- .tdr_cache_write(data$fuzzy.graphs[[sn]], cache_dir, "fuzzy.graphs", sn)
  }

  obj <- TDRObj(cellmap = cm)

  for (sn in samples) {
    expect_identical(
      .tdr_get_map_slot(obj, "clustering", sn), data$clustering[[sn]])
    expect_identical(
      .tdr_get_map_slot(obj, "celltyping", sn), data$celltyping[[sn]])
    expect_identical(
      .tdr_get_map_slot(obj, "nearest.lm", sn), data$nearest.lm[[sn]])
    expect_true(Matrix::all.equal(
      .tdr_get_map_slot(obj, "fuzzy.graphs", sn), data$fuzzy.graphs[[sn]]))
  }
})

# ── Backward compat: old slot names work ─────────────────────────────

test_that("legacy slot names map correctly via .tdr_get_map_slot", {
  obj <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = c(a = "1"))),
      celltyping   = list(ids = list(s1 = c(a = "A"))),
      nearest.lm   = list(s1 = matrix(1, 1, 1)),
      fuzzy.graphs = list(s1 = Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
    )
  )

  # Legacy names should all work
  expect_identical(.tdr_get_map_slot(obj, "clustering.ids", "s1"), c(a = "1"))
  expect_identical(.tdr_get_map_slot(obj, "celltyping.ids", "s1"), c(a = "A"))
  expect_identical(.tdr_get_map_slot(obj, "nearest.landmarks", "s1"), matrix(1, 1, 1))

  fg <- .tdr_get_map_slot(obj, "fuzzy.graph", "s1")
  expect_s4_class(fg, "sparseMatrix")
})

test_that("legacy slot names map correctly via .tdr_get_map_slot_all", {
  data_list <- list(s1 = c(a = "1"), s2 = c(b = "2"))
  obj <- TDRObj(
    cellmap = list(clustering = list(ids = data_list))
  )

  result <- .tdr_get_map_slot_all(obj, "clustering.ids")
  expect_identical(result, data_list)
})

# ── Validate / info / cleanup with new structure ─────────────────────

test_that("validate + info + cleanup full cycle with path strings", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_full_cycle")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering",   "s1")
  p_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graphs", "s1")

  obj <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  # Validate passes
  expect_message(tdr_cache_validate(obj), "Cache OK: 2")

  # Info reports correct stats
  expect_message(info <- tdr_cache_info(obj), "ACTIVE")
  expect_true(info$active)
  expect_equal(info$n_files, 2L)

  # Cleanup removes everything
  obj2 <- tdr_cache_cleanup(obj)
  expect_false(dir.exists(cache_dir))
  expect_null(obj2@cellmap$clustering$ids)
  expect_null(obj2@cellmap$fuzzy.graphs)
  expect_null(obj2@config$.cache.root)
})

# ── In-memory cellmap ignores validate ───────────────────────────────

test_that("in-memory cellmap entries are ignored by validate", {
  obj <- TDRObj(
    cellmap = list(
      clustering = list(ids = list(s1 = c(a = "cl.01"))),
      nearest.lm = list(s1 = matrix(1, 1, 1))
    )
  )

  # No path strings → "nothing to validate"
  expect_message(tdr_cache_validate(obj), "No on-disk cache active")
})

# ── Mixed in-memory / on-disk (not typical but should not crash) ─────

test_that("mixed in-memory and on-disk entries coexist safely", {
  cache_dir <- file.path(tempdir(), "tdr_test_cellmap_mixed")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  data_cl <- c(a = "cl.01")
  p_fg <- .tdr_cache_write("fg_data", cache_dir, "fuzzy.graphs", "s1")

  obj <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = data_cl)),   # in-memory
      fuzzy.graphs = list(s1 = p_fg)                    # on-disk
    )
  )

  # Both should be readable
  expect_identical(.tdr_get_map_slot(obj, "clustering", "s1"), data_cl)
  expect_identical(.tdr_get_map_slot(obj, "fuzzy.graphs", "s1"), "fg_data")
})

# ── get.cellmap() public API ────────────────────────────────────────

test_that("get.cellmap with new canonical slot names", {
  obj <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = c(a = "1", b = "2"))),
      celltyping   = list(ids = list(s1 = c(a = "A", b = "B"))),
      nearest.lm   = list(s1 = matrix(1:4, 2, 2)),
      fuzzy.graphs = list(s1 = Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
    )
  )

  expect_identical(get.cellmap(obj, "clustering", "s1"), c(a = "1", b = "2"))
  expect_identical(get.cellmap(obj, "celltyping", "s1"), c(a = "A", b = "B"))
  expect_identical(get.cellmap(obj, "nearest.lm", "s1"), matrix(1:4, 2, 2))
  expect_s4_class(get.cellmap(obj, "fuzzy.graphs", "s1"), "sparseMatrix")
})

test_that("get.cellmap backward-compat with legacy slot names", {
  obj <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = c(x = "cl1"))),
      celltyping   = list(ids = list(s1 = c(x = "ct1"))),
      nearest.lm   = list(s1 = matrix(9, 1, 1)),
      fuzzy.graphs = list(s1 = Matrix::sparseMatrix(i = 1, j = 1, x = 0.5, dims = c(1, 1)))
    )
  )

  # Legacy names should all work
  expect_identical(get.cellmap(obj, "clustering.ids", "s1"), c(x = "cl1"))
  expect_identical(get.cellmap(obj, "celltyping.ids", "s1"), c(x = "ct1"))
  expect_identical(get.cellmap(obj, "nearest.landmarks", "s1"), matrix(9, 1, 1))
  expect_s4_class(get.cellmap(obj, "fuzzy.graph", "s1"), "sparseMatrix")
})

test_that("get.cellmap errors on non-TDRObj input", {
  expect_error(get.cellmap(list(), "clustering", "s1"), "must be a TDRObj")
  expect_error(get.cellmap(NULL, "clustering", "s1"), "must be a TDRObj")
})

test_that("get.cellmap errors on non-existent sample", {
  obj <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = c(a = "1"))))
  )

  expect_error(
    get.cellmap(obj, "clustering", "no_such_sample"),
    "No entry for slot"
  )
})

# ── .annot_slot_map negative assertions ──────────────────────────────

test_that(".annot_slot_map returns only $slot and $comp_slot (no legacy keys)", {
  for (type in c("clustering", "celltyping")) {
    m <- .annot_slot_map(type)
    expect_named(m, c("slot", "comp_slot"), info = paste("type:", type))
    # Verify old keys are NOT present
    expect_null(m$cellmap_slot, info = paste("type:", type))
    expect_null(m$cache_slot, info = paste("type:", type))
  }
})

# ── Multi-sample validate/cleanup ───────────────────────────────────

test_that("validate + cleanup works with multiple samples", {
  cache_dir <- file.path(tempdir(), "tdr_test_multisample_cycle")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cm <- list(clustering = list(ids = list()), fuzzy.graphs = list())
  for (sn in c("s1", "s2", "s3")) {
    cm$clustering$ids[[sn]] <- .tdr_cache_write(
      c(a = paste0("cl.", sn)), cache_dir, "clustering", sn)
    cm$fuzzy.graphs[[sn]] <- .tdr_cache_write(
      Matrix::rsparsematrix(3, 2, 0.5), cache_dir, "fuzzy.graphs", sn)
  }

  obj <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = cm
  )

  # Validate passes for all 6 files
  expect_message(tdr_cache_validate(obj), "Cache OK: 6")

  # Info reports correct counts
  info <- capture.output(tdr_cache_info(obj), type = "message")
  expect_true(any(grepl("6", info)))

  # Cleanup removes all
  obj2 <- tdr_cache_cleanup(obj)
  expect_false(dir.exists(cache_dir))
  expect_null(obj2@cellmap$clustering$ids)
  expect_null(obj2@cellmap$fuzzy.graphs)
})

# ── Error paths ──────────────────────────────────────────────────────

test_that(".tdr_get_map_slot errors on invalid slot name", {
  obj <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = c(a = "1"))))
  )

  expect_error(
    .tdr_get_map_slot(obj, "invalid_slot", "s1"),
    "Unknown slot_name"
  )
})

test_that(".tdr_get_map_slot errors when cached file is missing", {
  fake_path <- file.path(tempdir(), "tdr_missing_file_test", "clustering", "s1.rds")
  attr(fake_path, "schema_v") <- .TDR_CACHE_SCHEMA_VERSION
  attr(fake_path, "bytes") <- 100L

  obj <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = fake_path)))
  )

  expect_error(
    .tdr_get_map_slot(obj, "clustering", "s1"),
    "Cache file missing"
  )
})
