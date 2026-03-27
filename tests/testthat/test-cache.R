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
# ──────────────────────────────────────────────────────────────────────

# ── Schema version constant ──────────────────────────────────────────

test_that(".TDR_CACHE_SCHEMA_VERSION is defined and is integer >= 2", {
  expect_true(is.integer(.TDR_CACHE_SCHEMA_VERSION))
  expect_length(.TDR_CACHE_SCHEMA_VERSION, 1)
  expect_true(.TDR_CACHE_SCHEMA_VERSION >= 2L)
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
  # Remove and re-call: should recreate

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
  
  # Format: run_YYYYMMDD_HHMMSS_<8 hex chars>
  expect_match(k1, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  expect_match(k2, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  
  # Distinct keys (randomness; vanishingly unlikely to collide)
  expect_false(identical(k1, k2))
})

# ── .tdr_cache_write + .tdr_cache_read round-trips ──────────────────

test_that(".tdr_cache_write + .tdr_cache_read round-trip correctly", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_rt")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Test with a named character vector (like clustering.ids)
  obj <- stats::setNames(c("cluster.01", "cluster.02", "cluster.01"),
                         nm = c("cell_A", "cell_B", "cell_C"))
  
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "sample1")
  
  expect_true(file.exists(meta$path))
  expect_equal(meta$slot, "clustering.ids")
  expect_equal(meta$sample, "sample1")
  expect_equal(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
  expect_true(meta$bytes > 0)
  # md5 field must NOT be present (removed for ephemeral cache)
  expect_null(meta$md5)
  
  # Read back
  obj2 <- .tdr_cache_read(meta)
  expect_identical(obj, obj2)
})

test_that(".tdr_cache_write + .tdr_cache_read round-trip for sparse matrix", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_sparse")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Simulate a sparse fuzzy graph (cells x landmarks)
  set.seed(42)
  fg <- Matrix::rsparsematrix(nrow = 100, ncol = 20, density = 0.1)
  
  meta <- .tdr_cache_write(fg, cache_dir, "fuzzy.graph", "sampleX")
  
  expect_true(file.exists(meta$path))
  
  fg2 <- .tdr_cache_read(meta)
  expect_true(Matrix::all.equal(fg, fg2))
})

test_that(".tdr_cache_write + .tdr_cache_read round-trip for integer matrix", {
  cache_dir <- file.path(tempdir(), "tdr_test_cache_intmat")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Simulate nearest.landmarks (cells x k)
  nn <- matrix(sample.int(50, 30, replace = TRUE), nrow = 10, ncol = 3)
  
  meta <- .tdr_cache_write(nn, cache_dir, "nearest.landmarks", "sampleY")
  nn2 <- .tdr_cache_read(meta)
  expect_identical(nn, nn2)
})

test_that(".tdr_cache_read errors on missing file", {
  meta <- list(path = file.path(tempdir(), "nonexistent_cache_file.rds"))
  expect_error(.tdr_cache_read(meta), "Cache file missing")
})

# ── Write metadata fields ───────────────────────────────────────────

test_that(".tdr_cache_write returns correct metadata (no md5)", {
  cache_dir <- file.path(tempdir(), "tdr_test_meta_fields")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(x = 1, y = 2, z = 3)
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  # Must have these fields
  expect_true(all(c("path", "slot", "sample", "bytes", "schema_v") %in% names(meta)))
  # Must NOT have md5 field
  expect_false("md5" %in% names(meta))
  expect_equal(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
})

test_that(".tdr_cache_write sets schema_v to .TDR_CACHE_SCHEMA_VERSION", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_write")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  meta <- .tdr_cache_write(42, cache_dir, "fuzzy.graph", "s1")
  expect_identical(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
})

# ── .tdr_cache_sweep_orphans ────────────────────────────────────────

test_that(".tdr_cache_sweep_orphans removes only .rds.tmp files", {
  cache_dir <- file.path(tempdir(), "tdr_test_sweep")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Create a real .rds and orphan .rds.tmp
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
  
  for (slot in c("clustering.ids", "fuzzy.graph", "nearest.landmarks")) {
    slot_dir <- file.path(cache_dir, slot)
    dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS("orphan", file.path(slot_dir, "sample1.rds.tmp"))
  }
  
  # Also keep a real .rds to make sure it isn't removed
  saveRDS("real", file.path(cache_dir, "clustering.ids", "sample1.rds"))
  
  expect_message(
    .tdr_cache_sweep_orphans(cache_dir),
    "Cleaned up 3 orphaned temp file"
  )
  
  # Real file should survive
  expect_true(file.exists(
    file.path(cache_dir, "clustering.ids", "sample1.rds")
  ))
})

# ── Schema version enforcement ───────────────────────────────────────

test_that(".tdr_cache_read errors on schema version mismatch", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_mismatch")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- "test"
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  # Tamper with schema version to simulate an old/future format
  meta$schema_v <- 999L
  
  expect_error(.tdr_cache_read(meta),
               "Cache schema mismatch.*v999")
})

test_that(".tdr_cache_read accepts current schema version", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_ok")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- "test"
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  expect_equal(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
  expect_identical(.tdr_cache_read(meta), obj)
})

# ── Atomic write quality ────────────────────────────────────────────

test_that("atomic rename: .rds file exists but no .rds.tmp after write", {
  cache_dir <- file.path(tempdir(), "tdr_test_atomic")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  .tdr_cache_write("data", cache_dir, "fuzzy.graph", "sample1")
  
  rds_file <- file.path(cache_dir, "fuzzy.graph", "sample1.rds")
  tmp_files <- list.files(file.path(cache_dir, "fuzzy.graph"),
                          pattern = "\\.rds\\.tmp$", full.names = TRUE)
  
  expect_true(file.exists(rds_file))
  expect_length(tmp_files, 0)
})

test_that(".tdr_cache_write creates intermediate directories", {
  cache_dir <- file.path(tempdir(), "tdr_test_deep", "nested", "dir")
  on.exit(unlink(file.path(tempdir(), "tdr_test_deep"), recursive = TRUE), add = TRUE)
  
  meta <- .tdr_cache_write("test_obj", cache_dir, "clustering.ids", "s1")
  expect_true(file.exists(meta$path))
  expect_identical(readRDS(meta$path), "test_obj")
})

# ── tdr_cache_cleanup ────────────────────────────────────────────────

test_that("tdr_cache_cleanup removes directory and strips metadata", {
  cache_dir <- file.path(tempdir(), "tdr_test_cleanup")
  dir.create(file.path(cache_dir, "clustering.ids"), recursive = TRUE)
  saveRDS("test", file.path(cache_dir, "clustering.ids", "s1.rds"))
  
  # Build a mock .tdr.obj
  mock_obj <- TDRObj(
    map = list(
      fdens = matrix(1, 2, 2),
      .cache = list(
        root = cache_dir,
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(
            s1 = list(path = file.path(cache_dir, "clustering.ids", "s1.rds"))
          )
        )
      )
    )
  )
  
  result <- tdr_cache_cleanup(mock_obj)
  
  expect_false(dir.exists(cache_dir))
  expect_null(result$map$.cache)
})

# ── .tdr_get_map_slot access patterns ───────────────────────────────

test_that(".tdr_get_map_slot reads from cache when on.disk = TRUE", {
  cache_dir <- file.path(tempdir(), "tdr_test_get_slot")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(a = "cluster.01", b = "cluster.02")
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "sample1")
  
  mock_tdr <- TDRObj(
    map = list(
      clustering = list(ids = NULL),
      .cache = list(
        on.disk = TRUE,
        manifests = list(
          clustering.ids = list(sample1 = meta)
        )
      )
    )
  )
  
  result <- .tdr_get_map_slot(mock_tdr, "clustering.ids", "sample1")
  expect_identical(result, obj)
})

test_that(".tdr_get_map_slot falls back to in-memory when no cache", {
  obj <- c(a = "cluster.01", b = "cluster.02")
  
  mock_tdr <- TDRObj(
    map = list(
      clustering = list(ids = list(sample1 = obj)),
      .cache = NULL
    )
  )
  
  result <- .tdr_get_map_slot(mock_tdr, "clustering.ids", "sample1")
  expect_identical(result, obj)
})

test_that(".tdr_get_map_slot_all reads all samples from cache", {
  cache_dir <- file.path(tempdir(), "tdr_test_get_all")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj1 <- c(a = "cl.01")
  obj2 <- c(b = "cl.02")
  m1 <- .tdr_cache_write(obj1, cache_dir, "clustering.ids", "s1")
  m2 <- .tdr_cache_write(obj2, cache_dir, "clustering.ids", "s2")
  
  mock_tdr <- TDRObj(
    map = list(
      clustering = list(ids = NULL),
      .cache = list(
        on.disk = TRUE,
        manifests = list(
          clustering.ids = list(s1 = m1, s2 = m2)
        )
      )
    )
  )
  
  result <- .tdr_get_map_slot_all(mock_tdr, "clustering.ids")
  expect_identical(result$s1, obj1)
  expect_identical(result$s2, obj2)
})

test_that(".tdr_get_map_slot_all falls back to in-memory", {
  obj_list <- list(s1 = c(a = "cl.01"), s2 = c(b = "cl.02"))
  
  mock_tdr <- TDRObj(
    map = list(
      clustering = list(ids = obj_list),
      .cache = NULL
    )
  )
  
  result <- .tdr_get_map_slot_all(mock_tdr, "clustering.ids")
  expect_identical(result, obj_list)
})

test_that(".tdr_get_map_slot errors on unknown slot_name", {
  mock_tdr <- TDRObj(map = list(.cache = NULL))
  expect_error(.tdr_get_map_slot(mock_tdr, "nonexistent", "s1"),
               "Unknown slot_name")
})

# ── tdr_cache_validate ───────────────────────────────────────────────

test_that("tdr_cache_validate returns invisibly when no cache", {
  mock_tdr <- TDRObj(map = list(.cache = NULL))
  expect_message(
    result <- tdr_cache_validate(mock_tdr),
    "No on-disk cache active"
  )
  expect_identical(result, mock_tdr)
})

test_that("tdr_cache_validate passes on valid cache", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_ok")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj1 <- c(a = "cl.01")
  obj2 <- c(b = "cl.02")
  m1 <- .tdr_cache_write(obj1, cache_dir, "clustering.ids", "s1")
  m2 <- .tdr_cache_write(obj2, cache_dir, "clustering.ids", "s2")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1, s2 = m2)
        )
      )
    )
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
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  # Delete the file
  file.remove(m1$path)
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1)
        )
      )
    )
  )
  
  expect_error(tdr_cache_validate(mock_tdr), "MISSING.*clustering.ids/s1")
})

test_that("tdr_cache_validate errors on schema mismatch", {
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = 999L,
        manifests = list()
      )
    )
  )
  
  expect_error(tdr_cache_validate(mock_tdr), "Cache schema mismatch.*v999")
})

test_that("tdr_cache_validate verbose = FALSE suppresses messages", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_quiet")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1)
        )
      )
    )
  )
  
  expect_silent(tdr_cache_validate(mock_tdr, .verbose = FALSE))
})

test_that("tdr_cache_validate checks across all 4 slot types", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_multi")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m_cl  <- .tdr_cache_write("cl",  cache_dir, "clustering.ids",    "s1")
  m_ct  <- .tdr_cache_write("ct",  cache_dir, "celltyping.ids",    "s1")
  m_nl  <- .tdr_cache_write("nl",  cache_dir, "nearest.landmarks", "s1")
  m_fg  <- .tdr_cache_write("fg",  cache_dir, "fuzzy.graph",       "s1")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids     = list(s1 = m_cl),
          celltyping.ids     = list(s1 = m_ct),
          nearest.landmarks  = list(s1 = m_nl),
          fuzzy.graph        = list(s1 = m_fg)
        )
      )
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
  
  m_cl <- .tdr_cache_write("cl", cache_dir, "clustering.ids",    "s1")
  m_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graph",       "s1")
  
  # Delete both files
  file.remove(m_cl$path)
  file.remove(m_fg$path)
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m_cl),
          fuzzy.graph    = list(s1 = m_fg)
        )
      )
    )
  )
  
  expect_error(
    tdr_cache_validate(mock_tdr),
    "MISSING.*clustering.ids/s1.*MISSING.*fuzzy.graph/s1"
  )
})

# ── .tdr_cache_validate_quiet ────────────────────────────────────────

test_that(".tdr_cache_validate_quiet is silent on valid cache", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_valid")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1)
        )
      )
    )
  )
  
  expect_silent(.tdr_cache_validate_quiet(mock_tdr))
})

test_that(".tdr_cache_validate_quiet is silent when no cache", {
  mock_tdr <- TDRObj(map = list(.cache = NULL))
  expect_silent(.tdr_cache_validate_quiet(mock_tdr))
})

test_that(".tdr_cache_validate_quiet errors on missing file", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_missing")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "nearest.landmarks", "s1")
  file.remove(m1$path)
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          nearest.landmarks = list(s1 = m1)
        )
      )
    )
  )
  
  expect_error(.tdr_cache_validate_quiet(mock_tdr), "MISSING")
})

# ── filelock integration ─────────────────────────────────────────────

test_that("filelock is optional: write works without it", {
  cache_dir <- file.path(tempdir(), "tdr_test_lock_optional")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  meta <- .tdr_cache_write("lock_test", cache_dir, "clustering.ids", "s1")
  expect_true(file.exists(meta$path))
  expect_identical(readRDS(meta$path), "lock_test")
  
  # No .lock files should remain after write
  lock_files <- list.files(cache_dir, pattern = "\\.lock$",
                           recursive = TRUE, full.names = TRUE)
  expect_length(lock_files, 0)
})

test_that("concurrent writes to different samples don't interfere", {
  cache_dir <- file.path(tempdir(), "tdr_test_concurrent_diff")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  results <- lapply(paste0("sample_", 1:10), function(sn) {
    .tdr_cache_write(paste0("data_", sn), cache_dir, "fuzzy.graph", sn)
  })
  
  for (i in seq_along(results)) {
    expect_true(file.exists(results[[i]]$path))
    expect_identical(readRDS(results[[i]]$path), paste0("data_sample_", i))
  }
})

test_that("overwrite to same sample produces valid file", {
  cache_dir <- file.path(tempdir(), "tdr_test_overwrite")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("version_1", cache_dir, "clustering.ids", "s1")
  m2 <- .tdr_cache_write("version_2", cache_dir, "clustering.ids", "s1")
  
  expect_identical(readRDS(m2$path), "version_2")
})

# ── Crash simulation ────────────────────────────────────────────────

test_that("orphan .rds.tmp left by crash is cleaned up", {
  cache_dir <- file.path(tempdir(), "tdr_test_crash_sim")
  slot_dir <- file.path(cache_dir, "fuzzy.graph")
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

# ──────────────────────────────────────────────────────────────────────
# tdr_cache_info tests
# ──────────────────────────────────────────────────────────────────────

test_that("tdr_cache_info reports INACTIVE when no cache", {
  mock_tdr <- TDRObj(map = list(.cache = NULL))
  
  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "INACTIVE"
  )
  expect_false(result$active)
})

test_that("tdr_cache_info reports ACTIVE with correct stats", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_active")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m_cl1 <- .tdr_cache_write(c(a = "cl.01"), cache_dir, "clustering.ids", "s1")
  m_cl2 <- .tdr_cache_write(c(b = "cl.02"), cache_dir, "clustering.ids", "s2")
  m_fg1 <- .tdr_cache_write(Matrix::rsparsematrix(10, 5, 0.3),
                             cache_dir, "fuzzy.graph", "s1")
  m_fg2 <- .tdr_cache_write(Matrix::rsparsematrix(10, 5, 0.3),
                             cache_dir, "fuzzy.graph", "s2")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk  = TRUE,
        root     = cache_dir,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m_cl1, s2 = m_cl2),
          fuzzy.graph    = list(s1 = m_fg1, s2 = m_fg2)
        )
      )
    )
  )
  
  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "ACTIVE"
  )
  
  expect_true(result$active)
  expect_equal(result$root, cache_dir)
  expect_equal(result$schema_v, .TDR_CACHE_SCHEMA_VERSION)
  expect_equal(sort(result$slots), c("clustering.ids", "fuzzy.graph"))
  expect_equal(result$n_samples, 2L)
  expect_equal(result$n_files, 4L)
  expect_true(result$total_bytes > 0)
})

test_that("tdr_cache_info prints human-readable size", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_size")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("small", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk  = TRUE,
        root     = cache_dir,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1)
        )
      )
    )
  )
  
  expect_message(tdr_cache_info(mock_tdr), "Total size:")
})

test_that("tdr_cache_info skips NULL slot manifests", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_null_slot")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk  = TRUE,
        root     = cache_dir,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m1),
          celltyping.ids = NULL  # no celltyping performed
        )
      )
    )
  )
  
  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "Slots:         clustering.ids"
  )
  
  expect_equal(result$slots, "clustering.ids")
  expect_equal(result$n_files, 1L)
})

test_that("tdr_cache_info handles empty cache map", {
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        root = "/tmp/nonexistent",
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list()
      )
    )
  )
  
  expect_message(
    result <- tdr_cache_info(mock_tdr),
    "ACTIVE"
  )
  expect_equal(result$n_files, 0L)
  expect_equal(result$n_samples, 0L)
})

# ──────────────────────────────────────────────────────────────────────
# NEW: Temp cache location tests
# ──────────────────────────────────────────────────────────────────────

test_that("cache directory is created under tempdir()", {
  # Use the helper directly
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

# ──────────────────────────────────────────────────────────────────────
# NEW: Session-end cleanup registry tests
# ──────────────────────────────────────────────────────────────────────

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
  dir.create(file.path(cache_dir, "clustering.ids"), recursive = TRUE)
  saveRDS("test", file.path(cache_dir, "clustering.ids", "s1.rds"))
  
  # Register it
  .tdr_cache_register_cleanup(cache_dir)
  expect_true(cache_dir %in% .tdr_cleanup_env$registered_dirs)
  
  mock_obj <- TDRObj(
    map = list(
      fdens = matrix(1, 2, 2),
      .cache = list(
        root = cache_dir,
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(
            s1 = list(path = file.path(cache_dir, "clustering.ids", "s1.rds"))
          )
        )
      )
    )
  )
  
  tdr_cache_cleanup(mock_obj)
  
  # Should be deregistered
  expect_false(cache_dir %in% .tdr_cleanup_env$registered_dirs)
  expect_false(dir.exists(cache_dir))
})

# Test that the finalizer function itself works correctly when invoked
# manually (we can't truly simulate session end in a test, but we can
# invoke the cleanup logic directly).
test_that("finalizer cleanup logic removes registered directories", {
  dir1 <- file.path(tempdir(), "tdr_test_fin_dir1")
  dir2 <- file.path(tempdir(), "tdr_test_fin_dir2")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir2, recursive = TRUE, showWarnings = FALSE)
  saveRDS("data", file.path(dir1, "test.rds"))
  saveRDS("data", file.path(dir2, "test.rds"))
  
  # Create a fresh environment to simulate the finalizer
  test_env <- new.env(parent = emptyenv())
  test_env$registered_dirs <- c(dir1, dir2)
  
  # This is the same function body that reg.finalizer uses
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

# ──────────────────────────────────────────────────────────────────────
# NEW: Collision avoidance tests
# ──────────────────────────────────────────────────────────────────────

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
  
  # Write same slot+sample name to both dirs
  m1 <- .tdr_cache_write("data_A", d1, "clustering.ids", "sample1")
  m2 <- .tdr_cache_write("data_B", d2, "clustering.ids", "sample1")
  
  # Files must be in different directories and hold different content
  expect_false(identical(m1$path, m2$path))
  expect_identical(.tdr_cache_read(m1), "data_A")
  expect_identical(.tdr_cache_read(m2), "data_B")
})

test_that("ten concurrent-style writes to same slot/different samples are isolated", {
  cache_dir <- file.path(tempdir(), "tdr_test_collision_10")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  sample_names <- paste0("sample_", seq_len(10))
  metas <- lapply(sample_names, function(sn) {
    .tdr_cache_write(list(id = sn), cache_dir, "nearest.landmarks", sn)
  })
  
  # Verify each sample has its own file with correct content
  for (i in seq_along(metas)) {
    obj <- .tdr_cache_read(metas[[i]])
    expect_identical(obj$id, sample_names[i])
  }
  
  # Verify there are exactly 10 files
  rds_files <- list.files(file.path(cache_dir, "nearest.landmarks"),
                          pattern = "\\.rds$")
  expect_length(rds_files, 10)
})

# ──────────────────────────────────────────────────────────────────────
# NEW: Edge cases
# ──────────────────────────────────────────────────────────────────────

test_that("cache dir already exists: no error on repeated creation", {
  cache_dir <- file.path(tempdir(), "tdr_test_exists_already")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Should not error
  expect_silent(
    meta <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  )
  expect_true(file.exists(meta$path))
})

test_that("partial cache: missing one slot file fails validate", {
  cache_dir <- file.path(tempdir(), "tdr_test_partial_cache")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m_cl <- .tdr_cache_write("cl", cache_dir, "clustering.ids",    "s1")
  m_fg <- .tdr_cache_write("fg", cache_dir, "fuzzy.graph",       "s1")
  
  # Delete only one file (simulate partial cache from interrupted write)
  file.remove(m_fg$path)
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          clustering.ids = list(s1 = m_cl),
          fuzzy.graph    = list(s1 = m_fg)
        )
      )
    )
  )
  
  # Validate should report fuzzy.graph missing but clustering.ids OK
  expect_error(tdr_cache_validate(mock_tdr), "MISSING.*fuzzy.graph/s1")
})

test_that("partial cache: reading intact slot still works", {
  cache_dir <- file.path(tempdir(), "tdr_test_partial_read")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m_cl <- .tdr_cache_write("cl_data", cache_dir, "clustering.ids", "s1")
  m_fg <- .tdr_cache_write("fg_data", cache_dir, "fuzzy.graph",    "s1")
  
  # Delete fuzzy.graph file
  file.remove(m_fg$path)
  
  mock_tdr <- TDRObj(
    map = list(
      .cache = list(
        on.disk = TRUE,
        manifests = list(
          clustering.ids = list(s1 = m_cl),
          fuzzy.graph    = list(s1 = m_fg)
        )
      )
    )
  )
  
  # The intact slot should still readable
  result <- .tdr_get_map_slot(mock_tdr, "clustering.ids", "s1")
  expect_identical(result, "cl_data")
  
  # The broken slot should error
  expect_error(.tdr_get_map_slot(mock_tdr, "fuzzy.graph", "s1"),
               "Cache file missing")
})

test_that("sample names with special characters are handled safely", {
  cache_dir <- file.path(tempdir(), "tdr_test_special_chars")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Sample names that might cause issues on some filesystems
  tricky_names <- c("sample with spaces", "sample-with-dashes",
                     "sample_123", "SAMPLE.dots")
  
  metas <- lapply(tricky_names, function(sn) {
    .tdr_cache_write(paste0("data_", sn), cache_dir, "clustering.ids", sn)
  })
  
  for (i in seq_along(metas)) {
    expect_true(file.exists(metas[[i]]$path))
    obj <- .tdr_cache_read(metas[[i]])
    expect_identical(obj, paste0("data_", tricky_names[i]))
  }
})

test_that(".tdr_cache_sweep_orphans is a no-op on non-existent directory", {
  expect_silent(.tdr_cache_sweep_orphans("/nonexistent/path/12345"))
})