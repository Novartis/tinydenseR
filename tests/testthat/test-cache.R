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
# Tests for on-disk caching infrastructure (R/cache.R)
#
# Cache is ephemeral — stored under tempdir() and cleaned up at
# session end via a registered finalizer.
#
# Schema v3: .tdr_cache_write returns attributed path strings;
# @cellmap stores path strings directly (no @density$.cache manifests).
# ──────────────────────────────────────────────────────────────────────

# ── Schema version constant ──────────────────────────────────────────

test_that(".TDR_CACHE_SCHEMA_VERSION is defined and is integer >= 3", {
  expect_true(is.integer(.TDR_CACHE_SCHEMA_VERSION))
  expect_length(.TDR_CACHE_SCHEMA_VERSION, 1)
  expect_true(.TDR_CACHE_SCHEMA_VERSION >= 3L)
})

# ── .tdr_cache_root ─────────────────────────────────────────────────

test_that(".tdr_cache_root returns a path under tempdir()", {
  root <- .tdr_cache_root()
  expect_true(startsWith(root, tempdir()))
  expect_true(grepl("tinydenseR_cache", root, fixed = TRUE))
  expect_true(dir.exists(root))
})

test_that(".tdr_cache_root creates directory lazily", {
  root <- .tdr_cache_root()
  unlink(root, recursive = TRUE)
  expect_false(dir.exists(root))
  root2 <- .tdr_cache_root()
  expect_true(dir.exists(root2))
  expect_identical(root, root2)
})

# ── .tdr_make_run_key ───────────────────────────────────────────────

test_that(".tdr_make_run_key generates unique, well-formatted keys", {
  k1 <- .tdr_make_run_key()
  k2 <- .tdr_make_run_key()
  expect_match(k1, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  expect_match(k2, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  expect_false(identical(k1, k2))
})

# ── .tdr_cache_write returns attributed path string ─────────────────

test_that(".tdr_cache_write returns attributed path string", {
  cache_dir <- file.path(tempdir(), "tdr_test_write_path")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  obj <- stats::setNames(c("cluster.01", "cluster.02", "cluster.01"),
                         nm = c("cell_A", "cell_B", "cell_C"))

  path <- .tdr_cache_write(obj, cache_dir, "clustering", "sample1")

  # Returns a character scalar

  expect_true(is.character(path))
  expect_length(path, 1)
  expect_true(file.exists(path))

  # Has schema_v and bytes attributes
  expect_identical(attr(path, "schema_v"), .TDR_CACHE_SCHEMA_VERSION)
  expect_true(attr(path, "bytes") > 0)

  # No list structure — just a string
  expect_false(is.list(path))
})

test_that(".tdr_cache_write + .tdr_cache_read round-trip correctly", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_rt")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  obj <- stats::setNames(c("cluster.01", "cluster.02", "cluster.01"),
                         nm = c("cell_A", "cell_B", "cell_C"))

  path <- .tdr_cache_write(obj, cache_dir, "clustering", "sample1")
  obj2 <- .tdr_cache_read(path)
  expect_identical(obj, obj2)
})

test_that(".tdr_cache_write + .tdr_cache_read round-trip for sparse matrix", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_sparse")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  set.seed(42)
  fg <- Matrix::rsparsematrix(nrow = 100, ncol = 20, density = 0.1)

  path <- .tdr_cache_write(fg, cache_dir, "fuzzy.graphs", "sampleX")
  expect_true(file.exists(path))

  fg2 <- .tdr_cache_read(path)
  expect_true(Matrix::all.equal(fg, fg2))
})

test_that(".tdr_cache_write + .tdr_cache_read round-trip for integer matrix", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_intmat")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  nn <- matrix(sample.int(50, 30, replace = TRUE), nrow = 10, ncol = 3)

  path <- .tdr_cache_write(nn, cache_dir, "nearest.lm", "sampleY")
  nn2 <- .tdr_cache_read(path)
  expect_identical(nn, nn2)
})

test_that(".tdr_cache_read errors on missing file", {
  path <- file.path(tempdir(), "nonexistent_cache_file.rds")
  attr(path, "schema_v") <- .TDR_CACHE_SCHEMA_VERSION
  expect_error(.tdr_cache_read(path), "Cache file missing")
})

test_that(".tdr_cache_read accepts legacy metadata list", {
  cache_dir <- file.path(tempdir(), "tdr_test_legacy_read")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  path <- .tdr_cache_write("test_obj", cache_dir, "clustering", "s1")

  # Build an old-style metadata list (schema v2 format)
  meta <- list(
    path = as.character(path),
    schema_v = .TDR_CACHE_SCHEMA_VERSION,
    bytes = 100
  )

  result <- .tdr_cache_read(meta)
  expect_identical(result, "test_obj")
})

# ── Path string attribute preservation ──────────────────────────────

test_that("path string attributes survive list operations", {
  cache_dir <- file.path(tempdir(), "tdr_test_attr_preserve")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  paths <- list()
  for (sn in c("s1", "s2", "s3")) {
    p <- .tdr_cache_write(paste0("data_", sn), cache_dir, "clustering", sn)
    paths[[sn]] <- p
  }

  # Verify attributes survive list storage
  for (sn in names(paths)) {
    expect_identical(attr(paths[[sn]], "schema_v"), .TDR_CACHE_SCHEMA_VERSION)
    expect_true(attr(paths[[sn]], "bytes") > 0)
  }
})

# ── Schema version enforcement ───────────────────────────────────────

test_that(".tdr_cache_read errors on schema version mismatch", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_mismatch")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  path <- .tdr_cache_write("test", cache_dir, "clustering", "s1")

  # Tamper with schema version
  attr(path, "schema_v") <- 999L
  expect_error(.tdr_cache_read(path), "Cache schema mismatch.*v999")
})

test_that(".tdr_cache_read accepts current schema version", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_ok")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  path <- .tdr_cache_write("test", cache_dir, "clustering", "s1")
  expect_identical(attr(path, "schema_v"), .TDR_CACHE_SCHEMA_VERSION)
  expect_identical(.tdr_cache_read(path), "test")
})

# ── .tdr_cache_sweep_orphans ────────────────────────────────────────

test_that(".tdr_cache_sweep_orphans removes only .rds.tmp files", {
  cache_dir <- file.path(tempdir(), "tdr_test_sweep")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  real <- file.path(cache_dir, "real.rds")
  orphan <- file.path(cache_dir, "orphan.rds.tmp")
  saveRDS("real", real)
  saveRDS("orphan", orphan)

  expect_true(file.exists(orphan))
  .tdr_cache_sweep_orphans(cache_dir)
  expect_false(file.exists(orphan))
  expect_true(file.exists(real))
})

test_that("sweep handles multiple orphans across slots", {
  cache_dir <- file.path(tempdir(), "tdr_test_multi_orphan")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  for (slot in c("clustering", "fuzzy.graphs", "nearest.lm")) {
    slot_dir <- file.path(cache_dir, slot)
    dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS("orphan", file.path(slot_dir, "sample1.rds.tmp"))
  }

  saveRDS("real", file.path(cache_dir, "clustering", "sample1.rds"))

  expect_message(
    .tdr_cache_sweep_orphans(cache_dir),
    "Cleaned up 3 orphaned temp file"
  )

  expect_true(file.exists(
    file.path(cache_dir, "clustering", "sample1.rds")
  ))
})

# ── Atomic write quality ────────────────────────────────────────────

test_that("atomic rename: .rds file exists but no .rds.tmp after write", {
  cache_dir <- file.path(tempdir(), "tdr_test_atomic")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  .tdr_cache_write("data", cache_dir, "fuzzy.graphs", "sample1")

  rds_file <- file.path(cache_dir, "fuzzy.graphs", "sample1.rds")
  tmp_files <- list.files(file.path(cache_dir, "fuzzy.graphs"),
                          pattern = "\\.rds\\.tmp$", full.names = TRUE)

  expect_true(file.exists(rds_file))
  expect_length(tmp_files, 0)
})

test_that(".tdr_cache_write creates intermediate directories", {
  cache_dir <- file.path(tempdir(), "tdr_test_deep", "nested", "dir")
  on.exit(unlink(file.path(tempdir(), "tdr_test_deep"), recursive = TRUE), add = TRUE)

  path <- .tdr_cache_write("test_obj", cache_dir, "clustering", "s1")
  expect_true(file.exists(path))
  expect_identical(readRDS(path), "test_obj")
})

# ── .tdr_get_map_slot with new cellmap structure ────────────────────

test_that(".tdr_get_map_slot reads in-memory data from new structure", {
  obj <- c(a = "cluster.01", b = "cluster.02")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(sample1 = obj)))
  )

  result <- .tdr_get_map_slot(mock_tdr, "clustering", "sample1")
  expect_identical(result, obj)
})

test_that(".tdr_get_map_slot auto-loads on-disk path strings", {
  cache_dir <- file.path(tempdir(), "tdr_test_get_slot_path")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  obj <- c(a = "cluster.01", b = "cluster.02")
  path <- .tdr_cache_write(obj, cache_dir, "clustering", "sample1")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(sample1 = path)))
  )

  result <- .tdr_get_map_slot(mock_tdr, "clustering", "sample1")
  expect_identical(result, obj)
})

test_that(".tdr_get_map_slot accepts legacy slot names", {
  obj <- c(a = "cluster.01")
  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = obj)))
  )

  # Legacy name "clustering.ids" → "clustering"
  result <- .tdr_get_map_slot(mock_tdr, "clustering.ids", "s1")
  expect_identical(result, obj)
})

test_that(".tdr_get_map_slot_all reads all samples", {
  cache_dir <- file.path(tempdir(), "tdr_test_get_all")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  obj1 <- c(a = "cl.01")
  obj2 <- c(b = "cl.02")
  p1 <- .tdr_cache_write(obj1, cache_dir, "clustering", "s1")
  p2 <- .tdr_cache_write(obj2, cache_dir, "clustering", "s2")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1, s2 = p2)))
  )

  result <- .tdr_get_map_slot_all(mock_tdr, "clustering")
  expect_identical(result$s1, obj1)
  expect_identical(result$s2, obj2)
})

test_that(".tdr_get_map_slot_all falls back to in-memory", {
  obj_list <- list(s1 = c(a = "cl.01"), s2 = c(b = "cl.02"))

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = obj_list))
  )

  result <- .tdr_get_map_slot_all(mock_tdr, "clustering")
  expect_identical(result, obj_list)
})

test_that(".tdr_get_map_slot errors on unknown slot_name", {
  mock_tdr <- TDRObj()
  expect_error(.tdr_get_map_slot(mock_tdr, "nonexistent", "s1"),
               "Unknown slot_name")
})

# ── tdr_cache_validate with path strings ─────────────────────────────

test_that("tdr_cache_validate returns invisibly when no cache", {
  mock_tdr <- TDRObj()
  expect_message(
    result <- tdr_cache_validate(mock_tdr),
    "No on-disk cache active"
  )
  expect_identical(result, mock_tdr)
})

test_that("tdr_cache_validate passes on valid path strings", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_ok")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("cl1", cache_dir, "clustering", "s1")
  p2 <- .tdr_cache_write("cl2", cache_dir, "clustering", "s2")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1, s2 = p2)))
  )

  expect_message(
    result <- tdr_cache_validate(mock_tdr),
    "Cache OK: 2 file\\(s\\) verified"
  )
  expect_identical(result, mock_tdr)
})

test_that("tdr_cache_validate errors on missing file", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_missing")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "clustering", "s1")
  file.remove(p1)

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  expect_error(tdr_cache_validate(mock_tdr), "MISSING.*clustering/s1")
})

test_that("tdr_cache_validate errors on schema mismatch", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_schema")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "clustering", "s1")
  attr(p1, "schema_v") <- 999L

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  expect_error(tdr_cache_validate(mock_tdr), "Cache schema mismatch.*v999")
})

test_that("tdr_cache_validate verbose = FALSE suppresses messages", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_quiet")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "clustering", "s1")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  expect_silent(tdr_cache_validate(mock_tdr, .verbose = FALSE))
})

test_that("tdr_cache_validate checks across all 4 slot types", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_multi")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering",   "s1")
  p_ct <- .tdr_cache_write("ct", cache_dir, "celltyping",   "s1")
  p_nl <- .tdr_cache_write("nl", cache_dir, "nearest.lm",   "s1")
  p_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graphs", "s1")

  mock_tdr <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      celltyping   = list(ids = list(s1 = p_ct)),
      nearest.lm   = list(s1 = p_nl),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  expect_message(
    tdr_cache_validate(mock_tdr),
    "Cache OK: 4 file\\(s\\) verified"
  )
})

test_that("tdr_cache_validate reports multiple missing files", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_multi_miss")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering",   "s1")
  p_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graphs", "s1")

  file.remove(p_cl)
  file.remove(p_fg)

  mock_tdr <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  expect_error(
    tdr_cache_validate(mock_tdr),
    "MISSING.*clustering/s1.*MISSING.*fuzzy.graphs/s1"
  )
})

# ── .tdr_cache_validate_quiet ────────────────────────────────────────

test_that(".tdr_cache_validate_quiet is silent on valid cache", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_valid")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "clustering", "s1")

  mock_tdr <- TDRObj(
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  expect_silent(.tdr_cache_validate_quiet(mock_tdr))
})

test_that(".tdr_cache_validate_quiet is silent when no cache", {
  mock_tdr <- TDRObj()
  expect_silent(.tdr_cache_validate_quiet(mock_tdr))
})

test_that(".tdr_cache_validate_quiet errors on missing file", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_missing")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "nearest.lm", "s1")
  file.remove(p1)

  mock_tdr <- TDRObj(
    cellmap = list(nearest.lm = list(s1 = p1))
  )

  expect_error(.tdr_cache_validate_quiet(mock_tdr), "MISSING")
})

# ── tdr_cache_cleanup ────────────────────────────────────────────────

test_that("tdr_cache_cleanup removes directory and clears cellmap", {
  cache_dir <- file.path(tempdir(), "tdr_test_cleanup")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering", "s1")

  mock_obj <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(clustering = list(ids = list(s1 = p_cl)))
  )

  result <- tdr_cache_cleanup(mock_obj)

  expect_false(dir.exists(cache_dir))
  expect_null(result@cellmap$clustering$ids)
  expect_null(result@config$.cache.root)
})

test_that("tdr_cache_cleanup derives root from path strings as fallback", {
  cache_dir <- file.path(tempdir(), "tdr_test_cleanup_fallback")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering", "s1")

  # No .cache.root in config
  mock_obj <- TDRObj(
    config = list(assay.type = "RNA"),
    cellmap = list(clustering = list(ids = list(s1 = p_cl)))
  )

  result <- tdr_cache_cleanup(mock_obj)
  expect_false(dir.exists(cache_dir))
})

# ── tdr_cache_info ───────────────────────────────────────────────────

test_that("tdr_cache_info reports INACTIVE when no cache", {
  mock_tdr <- TDRObj()

  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "INACTIVE"
  )
  expect_false(result$active)
})

test_that("tdr_cache_info reports ACTIVE with correct stats", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_active")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl1 <- .tdr_cache_write(c(a = "cl.01"), cache_dir, "clustering", "s1")
  p_cl2 <- .tdr_cache_write(c(b = "cl.02"), cache_dir, "clustering", "s2")
  p_fg1 <- .tdr_cache_write(Matrix::rsparsematrix(10, 5, 0.3),
                             cache_dir, "fuzzy.graphs", "s1")
  p_fg2 <- .tdr_cache_write(Matrix::rsparsematrix(10, 5, 0.3),
                             cache_dir, "fuzzy.graphs", "s2")

  mock_tdr <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl1, s2 = p_cl2)),
      fuzzy.graphs = list(s1 = p_fg1, s2 = p_fg2)
    )
  )

  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "ACTIVE"
  )

  expect_true(result$active)
  expect_equal(result$root, cache_dir)
  expect_equal(result$schema_v, .TDR_CACHE_SCHEMA_VERSION)
  expect_equal(sort(result$slots), c("clustering", "fuzzy.graphs"))
  expect_equal(result$n_samples, 2L)
  expect_equal(result$n_files, 4L)
  expect_true(result$total_bytes > 0)
})

test_that("tdr_cache_info prints human-readable size", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_size")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("small", cache_dir, "clustering", "s1")

  mock_tdr <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(clustering = list(ids = list(s1 = p1)))
  )

  expect_message(tdr_cache_info(mock_tdr), "Total size:")
})

test_that("tdr_cache_info skips NULL slot entries", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_null_slot")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("data", cache_dir, "clustering", "s1")

  mock_tdr <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(
      clustering = list(ids = list(s1 = p1)),
      celltyping = NULL
    )
  )

  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "Slots:         clustering"
  )

  expect_equal(result$slots, "clustering")
  expect_equal(result$n_files, 1L)
})

# ── filelock integration ─────────────────────────────────────────────

test_that("filelock is optional: write works without it", {
  cache_dir <- file.path(tempdir(), "tdr_test_lock_optional")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  path <- .tdr_cache_write("lock_test", cache_dir, "clustering", "s1")
  expect_true(file.exists(path))
  expect_identical(readRDS(path), "lock_test")

  lock_files <- list.files(cache_dir, pattern = "\\.lock$",
                           recursive = TRUE, full.names = TRUE)
  expect_length(lock_files, 0)
})

test_that("concurrent writes to different samples don't interfere", {
  cache_dir <- file.path(tempdir(), "tdr_test_concurrent_diff")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  results <- lapply(paste0("sample_", 1:10), function(sn) {
    .tdr_cache_write(paste0("data_", sn), cache_dir, "fuzzy.graphs", sn)
  })

  for (i in seq_along(results)) {
    expect_true(file.exists(results[[i]]))
    expect_identical(readRDS(results[[i]]), paste0("data_sample_", i))
  }
})

test_that("overwrite to same sample produces valid file", {
  cache_dir <- file.path(tempdir(), "tdr_test_overwrite")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p1 <- .tdr_cache_write("version_1", cache_dir, "clustering", "s1")
  p2 <- .tdr_cache_write("version_2", cache_dir, "clustering", "s1")

  expect_identical(readRDS(p2), "version_2")
})

# ── Session-end cleanup registry ─────────────────────────────────────

test_that(".tdr_cache_register_cleanup adds dir to registry", {
  test_dir <- file.path(tempdir(), "tdr_test_register_cleanup")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    unlink(test_dir, recursive = TRUE)
    .tdr_cache_deregister(test_dir)
  }, add = TRUE)

  .tdr_cache_register_cleanup(test_dir)
  expect_true(test_dir %in% .tdr_cleanup_env$registered_dirs)
})

test_that(".tdr_cache_deregister removes dir from registry", {
  test_dir <- file.path(tempdir(), "tdr_test_deregister")
  .tdr_cache_register_cleanup(test_dir)
  expect_true(test_dir %in% .tdr_cleanup_env$registered_dirs)

  .tdr_cache_deregister(test_dir)
  expect_false(test_dir %in% .tdr_cleanup_env$registered_dirs)
})

test_that("cleanup finalizer is marked as set after first registration", {
  test_dir <- file.path(tempdir(), "tdr_test_finalizer_flag")
  on.exit(.tdr_cache_deregister(test_dir), add = TRUE)

  .tdr_cache_register_cleanup(test_dir)
  expect_true(.tdr_cleanup_env$finalizer_set)
})

test_that("duplicate registration does not create duplicate entries", {
  test_dir <- file.path(tempdir(), "tdr_test_dup_register")
  on.exit(.tdr_cache_deregister(test_dir), add = TRUE)

  .tdr_cache_register_cleanup(test_dir)
  .tdr_cache_register_cleanup(test_dir)

  n_matches <- sum(.tdr_cleanup_env$registered_dirs == test_dir)
  expect_equal(n_matches, 1L)
})

test_that("tdr_cache_cleanup deregisters the directory", {
  cache_dir <- file.path(tempdir(), "tdr_test_cleanup_dereg")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("test", cache_dir, "clustering", "s1")
  .tdr_cache_register_cleanup(cache_dir)
  expect_true(cache_dir %in% .tdr_cleanup_env$registered_dirs)

  mock_obj <- TDRObj(
    config = list(assay.type = "RNA", .cache.root = cache_dir),
    cellmap = list(clustering = list(ids = list(s1 = p_cl)))
  )

  tdr_cache_cleanup(mock_obj)
  expect_false(cache_dir %in% .tdr_cleanup_env$registered_dirs)
  expect_false(dir.exists(cache_dir))
})

test_that("finalizer cleanup logic removes registered directories", {
  dir1 <- file.path(tempdir(), "tdr_test_fin_dir1")
  dir2 <- file.path(tempdir(), "tdr_test_fin_dir2")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir2, recursive = TRUE, showWarnings = FALSE)
  saveRDS("data", file.path(dir1, "test.rds"))
  saveRDS("data", file.path(dir2, "test.rds"))

  test_env <- new.env(parent = emptyenv())
  test_env$registered_dirs <- c(dir1, dir2)

  cleanup_fn <- function(env) {
    dirs <- env$registered_dirs
    for (d in dirs) {
      if (dir.exists(d)) {
        unlink(d, recursive = TRUE, force = TRUE)
      }
    }
  }

  cleanup_fn(test_env)
  expect_false(dir.exists(dir1))
  expect_false(dir.exists(dir2))
})

# ── Collision avoidance ──────────────────────────────────────────────

test_that("two independent run keys produce separate cache dirs", {
  root <- .tdr_cache_root()
  k1 <- .tdr_make_run_key()
  k2 <- .tdr_make_run_key()

  d1 <- file.path(root, k1)
  d2 <- file.path(root, k2)
  dir.create(d1, recursive = TRUE, showWarnings = FALSE)
  dir.create(d2, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    unlink(d1, recursive = TRUE)
    unlink(d2, recursive = TRUE)
  }, add = TRUE)

  p1 <- .tdr_cache_write("data_A", d1, "clustering", "sample1")
  p2 <- .tdr_cache_write("data_B", d2, "clustering", "sample1")

  expect_false(identical(p1, p2))
  expect_identical(.tdr_cache_read(p1), "data_A")
  expect_identical(.tdr_cache_read(p2), "data_B")
})

test_that("ten concurrent-style writes to same slot/different samples are isolated", {
  cache_dir <- file.path(tempdir(), "tdr_test_collision_10")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  sample_names <- paste0("sample_", seq_len(10))
  paths <- lapply(sample_names, function(sn) {
    .tdr_cache_write(list(id = sn), cache_dir, "nearest.lm", sn)
  })

  for (i in seq_along(paths)) {
    obj <- .tdr_cache_read(paths[[i]])
    expect_identical(obj$id, sample_names[i])
  }

  rds_files <- list.files(file.path(cache_dir, "nearest.lm"),
                          pattern = "\\.rds$")
  expect_length(rds_files, 10)
})

# ── Edge cases ───────────────────────────────────────────────────────

test_that("cache dir already exists: no error on repeated creation", {
  cache_dir <- file.path(tempdir(), "tdr_test_exists_already")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  expect_silent(
    path <- .tdr_cache_write("data", cache_dir, "clustering", "s1")
  )
  expect_true(file.exists(path))
})

test_that("partial cache: missing one slot file fails validate", {
  cache_dir <- file.path(tempdir(), "tdr_test_partial_cache")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl", cache_dir, "clustering",   "s1")
  p_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graphs", "s1")

  file.remove(p_fg)

  mock_tdr <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  expect_error(tdr_cache_validate(mock_tdr), "MISSING.*fuzzy.graphs/s1")
})

test_that("partial cache: reading intact slot still works", {
  cache_dir <- file.path(tempdir(), "tdr_test_partial_read")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  p_cl <- .tdr_cache_write("cl_data", cache_dir, "clustering",   "s1")
  p_fg <- .tdr_cache_write("fg_data", cache_dir, "fuzzy.graphs", "s1")

  file.remove(p_fg)

  mock_tdr <- TDRObj(
    cellmap = list(
      clustering   = list(ids = list(s1 = p_cl)),
      fuzzy.graphs = list(s1 = p_fg)
    )
  )

  result <- .tdr_get_map_slot(mock_tdr, "clustering", "s1")
  expect_identical(result, "cl_data")

  expect_error(.tdr_get_map_slot(mock_tdr, "fuzzy.graphs", "s1"),
               "Cache file missing")
})

test_that("sample names with special characters are handled safely", {
  cache_dir <- file.path(tempdir(), "tdr_test_special_chars")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  tricky_names <- c("sample with spaces", "sample-with-dashes",
                     "sample_123", "SAMPLE.dots")

  paths <- lapply(tricky_names, function(sn) {
    .tdr_cache_write(paste0("data_", sn), cache_dir, "clustering", sn)
  })

  for (i in seq_along(paths)) {
    expect_true(file.exists(paths[[i]]))
    obj <- .tdr_cache_read(paths[[i]])
    expect_identical(obj, paste0("data_", tricky_names[i]))
  }
})

test_that(".tdr_cache_sweep_orphans is a no-op on non-existent directory", {
  expect_silent(.tdr_cache_sweep_orphans("/nonexistent/path/12345"))
})

# ── Crash simulation ────────────────────────────────────────────────

test_that("orphan .rds.tmp left by crash is cleaned up", {
  cache_dir <- file.path(tempdir(), "tdr_test_crash_sim")
  slot_dir <- file.path(cache_dir, "fuzzy.graphs")
  dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  orphan <- file.path(slot_dir, "sample1.rds.tmp")
  saveRDS("crash_data", orphan)
  expect_true(file.exists(orphan))

  expect_message(
    .tdr_cache_sweep_orphans(cache_dir),
    "Cleaned up 1 orphaned temp file"
  )
  expect_false(file.exists(orphan))
})

# ── Backward compat: constructor migrates old @density$.cache ────────

test_that("TDRObj constructor migrates old @density$.cache to @cellmap path strings", {
  cache_dir <- file.path(tempdir(), "tdr_test_migration")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  # Write actual files
  dir.create(file.path(cache_dir, "clustering.ids"), recursive = TRUE)
  dir.create(file.path(cache_dir, "nearest.landmarks"), recursive = TRUE)
  cl_path <- file.path(cache_dir, "clustering.ids", "s1.rds")
  nl_path <- file.path(cache_dir, "nearest.landmarks", "s1.rds")
  saveRDS("cl_data", cl_path)
  saveRDS("nl_data", nl_path)

  # Build with old-style map containing .cache
  obj <- TDRObj(
    config = list(assay.type = "RNA"),
    map = list(
      norm = matrix(1, 2, 2),
      .cache = list(
        root = cache_dir,
        on.disk = TRUE,
        schema_v = 3L,
        manifests = list(
          clustering.ids = list(
            s1 = list(path = cl_path, schema_v = 3L, bytes = file.size(cl_path))
          ),
          nearest.landmarks = list(
            s1 = list(path = nl_path, schema_v = 3L, bytes = file.size(nl_path))
          )
        )
      )
    )
  )

  # @density$.cache should be gone
  expect_null(obj@density$.cache)

  # Path strings should be in @cellmap
  expect_true(!is.null(obj@cellmap$clustering$ids$s1))
  expect_true(!is.null(obj@cellmap$nearest.lm$s1))

  # Should be attributed path strings
  expect_identical(attr(obj@cellmap$clustering$ids$s1, "schema_v"), 3L)

  # Cache root stored in config
  expect_equal(obj@config$.cache.root, cache_dir)
})

test_that("TDRObj constructor migrates old cellmap slot names", {
  old_ids <- list(s1 = c(a = "cl.01"), s2 = c(b = "cl.02"))

  obj <- TDRObj(
    cellmap = list(
      cluster.ids  = old_ids,
      celltype.ids = old_ids,
      fuzzy.graph  = list(s1 = "fg1", s2 = "fg2"),
      nearest.lm   = list(s1 = "nl1", s2 = "nl2")
    )
  )

  # Old names should be migrated
  expect_identical(obj@cellmap$clustering$ids, old_ids)
  expect_identical(obj@cellmap$celltyping$ids, old_ids)
  expect_identical(obj@cellmap$fuzzy.graphs, list(s1 = "fg1", s2 = "fg2"))
  # nearest.lm stays as-is (no rename needed)
  expect_identical(obj@cellmap$nearest.lm, list(s1 = "nl1", s2 = "nl2"))

  # Old names should be cleaned up
  expect_null(obj@cellmap[["cluster.ids"]])
  expect_null(obj@cellmap[["celltype.ids"]])
  expect_null(obj@cellmap[["fuzzy.graph"]])
})

# ── Temp cache location tests ───────────────────────────────────────

test_that("cache directory is created under tempdir()", {
  root <- .tdr_cache_root()
  expect_true(startsWith(normalizePath(root, mustWork = FALSE),
                         normalizePath(tempdir(), mustWork = FALSE)))
})

test_that("run cache dirs are nested under the session cache root", {
  root <- .tdr_cache_root()
  run_key <- .tdr_make_run_key()
  run_dir <- file.path(root, run_key)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  expect_true(dir.exists(run_dir))
  expect_true(startsWith(run_dir, root))
})
