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
# ──────────────────────────────────────────────────────────────────────

test_that(".tdr_make_run_key generates unique, well-formatted keys", {
  k1 <- .tdr_make_run_key()
  k2 <- .tdr_make_run_key()
  
  # Format: run_YYYYMMDD_HHMMSS_<8 hex chars>

  expect_match(k1, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  expect_match(k2, "^run_\\d{8}_\\d{6}_[0-9a-f]{8}$")
  
  # Distinct keys (randomness; vanishingly unlikely to collide)
  expect_false(identical(k1, k2))
})

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
  expect_equal(meta$schema_v, 1L)
  expect_true(meta$bytes > 0)
  
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

test_that("tdr_cache_cleanup removes directory and strips metadata", {
  cache_dir <- file.path(tempdir(), "tdr_test_cleanup")
  dir.create(file.path(cache_dir, "clustering.ids"), recursive = TRUE)
  saveRDS("test", file.path(cache_dir, "clustering.ids", "s1.rds"))
  
  # Build a mock .tdr.obj
  mock_obj <- list(
    map = list(
      fdens = matrix(1, 2, 2),
      .cache = list(
        root = cache_dir,
        on.disk = TRUE,
        schema_v = 1L,
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

test_that(".tdr_get_map_slot reads from cache when on.disk = TRUE", {
  cache_dir <- file.path(tempdir(), "tdr_test_get_slot")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(a = "cluster.01", b = "cluster.02")
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "sample1")
  
  mock_tdr <- list(
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
  
  mock_tdr <- list(
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
  
  mock_tdr <- list(
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
  
  mock_tdr <- list(
    map = list(
      clustering = list(ids = obj_list),
      .cache = NULL
    )
  )
  
  result <- .tdr_get_map_slot_all(mock_tdr, "clustering.ids")
  expect_identical(result, obj_list)
})

test_that(".tdr_get_map_slot errors on unknown slot_name", {
  mock_tdr <- list(map = list(.cache = NULL))
  expect_error(.tdr_get_map_slot(mock_tdr, "nonexistent", "s1"),
               "Unknown slot_name")
})

test_that(".tdr_cache_write creates intermediate directories", {
  cache_dir <- file.path(tempdir(), "tdr_test_deep", "nested", "dir")
  on.exit(unlink(file.path(tempdir(), "tdr_test_deep"), recursive = TRUE), add = TRUE)
  
  meta <- .tdr_cache_write("test_obj", cache_dir, "clustering.ids", "s1")
  expect_true(file.exists(meta$path))
  expect_identical(readRDS(meta$path), "test_obj")
})

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


# ──────────────────────────────────────────────────────────────────────
# Phase 2 — Hardening tests
# ──────────────────────────────────────────────────────────────────────

# ── Schema version constant ──────────────────────────────────────────

test_that(".TDR_CACHE_SCHEMA_VERSION is defined and is integer", {
  expect_true(is.integer(.TDR_CACHE_SCHEMA_VERSION))
  expect_length(.TDR_CACHE_SCHEMA_VERSION, 1)
  expect_true(.TDR_CACHE_SCHEMA_VERSION >= 1L)
})

# ── MD5 checksum round-trip ──────────────────────────────────────────

test_that(".tdr_cache_write stores md5 in metadata", {
  cache_dir <- file.path(tempdir(), "tdr_test_md5_meta")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(x = 1, y = 2, z = 3)
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  expect_true(!is.null(meta$md5))
  expect_match(meta$md5, "^[0-9a-f]{32}$")
  
  # md5 should match the file on disk
  expect_equal(meta$md5, unname(tools::md5sum(meta$path)))
})

test_that(".tdr_cache_read with verify = TRUE succeeds on clean file", {
  cache_dir <- file.path(tempdir(), "tdr_test_md5_ok")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- list(a = 1:5, b = "hello")
  meta <- .tdr_cache_write(obj, cache_dir, "fuzzy.graph", "sA")
  
  result <- .tdr_cache_read(meta, verify = TRUE)
  expect_identical(result, obj)
})

test_that(".tdr_cache_read with verify = TRUE catches corruption", {
  cache_dir <- file.path(tempdir(), "tdr_test_md5_bad")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(a = 1, b = 2)
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  # Corrupt the file by overwriting with different content
  saveRDS("corrupted!", meta$path)
  
  expect_error(.tdr_cache_read(meta, verify = TRUE),
               "Cache checksum mismatch")
})

test_that(".tdr_cache_read without verify ignores checksum", {
  cache_dir <- file.path(tempdir(), "tdr_test_md5_noverify")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  obj <- c(a = 1, b = 2)
  meta <- .tdr_cache_write(obj, cache_dir, "clustering.ids", "s1")
  
  # Corrupt the file
  saveRDS("corrupted!", meta$path)
  
  # Should not error (verify = FALSE by default)
  result <- .tdr_cache_read(meta)
  expect_identical(result, "corrupted!")
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
  
  # schema_v should equal the constant

  expect_equal(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
  expect_identical(.tdr_cache_read(meta), obj)
})

test_that(".tdr_cache_write sets schema_v to .TDR_CACHE_SCHEMA_VERSION", {
  cache_dir <- file.path(tempdir(), "tdr_test_schema_write")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  meta <- .tdr_cache_write(42, cache_dir, "fuzzy.graph", "s1")
  expect_identical(meta$schema_v, .TDR_CACHE_SCHEMA_VERSION)
})

# ── tdr_cache_validate ───────────────────────────────────────────────

test_that("tdr_cache_validate returns invisibly when no cache", {
  mock_tdr <- list(map = list(.cache = NULL))
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
  
  mock_tdr <- list(
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
  
  mock_tdr <- list(
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

test_that("tdr_cache_validate with .verify.checksum catches corruption", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_checksum")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("good_data", cache_dir, "fuzzy.graph", "s1")
  
  # Corrupt the file
  saveRDS("bad_data", m1$path)
  
  mock_tdr <- list(
    map = list(
      .cache = list(
        on.disk = TRUE,
        schema_v = .TDR_CACHE_SCHEMA_VERSION,
        manifests = list(
          fuzzy.graph = list(s1 = m1)
        )
      )
    )
  )
  
  # Without checksum — should pass (file exists)
  expect_message(
    tdr_cache_validate(mock_tdr, .verify.checksum = FALSE),
    "Cache OK"
  )
  
  # With checksum — should fail
  expect_error(
    tdr_cache_validate(mock_tdr, .verify.checksum = TRUE),
    "CHECKSUM MISMATCH"
  )
})

test_that("tdr_cache_validate errors on schema mismatch", {
  mock_tdr <- list(
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
  
  mock_tdr <- list(
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

# ── .tdr_cache_validate_quiet ────────────────────────────────────────

test_that(".tdr_cache_validate_quiet is silent on valid cache", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_valid")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- list(
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
  mock_tdr <- list(map = list(.cache = NULL))
  expect_silent(.tdr_cache_validate_quiet(mock_tdr))
})

test_that(".tdr_cache_validate_quiet errors on missing file", {
  cache_dir <- file.path(tempdir(), "tdr_test_quiet_missing")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "nearest.landmarks", "s1")
  file.remove(m1$path)
  
  mock_tdr <- list(
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
  # This test always passes since filelock is Suggests, not Imports.
  # If filelock IS installed, the write uses it; if not, it skips locking.
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
  
  # Simulate sequential (not truly parallel, but checks no file clobbering)
  results <- lapply(paste0("sample_", 1:10), function(sn) {
    .tdr_cache_write(paste0("data_", sn), cache_dir, "fuzzy.graph", sn)
  })
  
  # All files should exist with correct content
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
  
  # Should have latest content
  expect_identical(readRDS(m2$path), "version_2")
  
  # md5 should differ
  expect_false(identical(m1$md5, m2$md5))
})

# ── Crash simulation ────────────────────────────────────────────────

test_that("orphan .rds.tmp left by crash is cleaned up", {
  cache_dir <- file.path(tempdir(), "tdr_test_crash_sim")
  slot_dir <- file.path(cache_dir, "fuzzy.graph")
  dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  # Simulate a crash: .rds.tmp exists but no .rds
  orphan <- file.path(slot_dir, "sample1.rds.tmp")
  saveRDS("crash_data", orphan)
  expect_true(file.exists(orphan))
  
  expect_message(
    .tdr_cache_sweep_orphans(cache_dir),
    "Cleaned up 1 orphaned temp file"
  )
  expect_false(file.exists(orphan))
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

# ── Multi-slot validate ──────────────────────────────────────────────

test_that("tdr_cache_validate checks across all 4 slot types", {
  cache_dir <- file.path(tempdir(), "tdr_test_validate_multi")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m_cl  <- .tdr_cache_write("cl",  cache_dir, "clustering.ids",    "s1")
  m_ct  <- .tdr_cache_write("ct",  cache_dir, "celltyping.ids",    "s1")
  m_nl  <- .tdr_cache_write("nl",  cache_dir, "nearest.landmarks", "s1")
  m_fg  <- .tdr_cache_write("fg",  cache_dir, "fuzzy.graph",       "s1")
  
  mock_tdr <- list(
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
  
  # With checksum
  expect_message(
    tdr_cache_validate(mock_tdr, .verify.checksum = TRUE),
    "with MD5"
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
  
  mock_tdr <- list(
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


# ──────────────────────────────────────────────────────────────────────
# Phase 3 — tdr_cache_info tests
# ──────────────────────────────────────────────────────────────────────

test_that("tdr_cache_info reports INACTIVE when no cache", {
  mock_tdr <- list(map = list(.cache = NULL))
  
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
  
  mock_tdr <- list(
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
  
  mock_tdr <- list(
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
  
  # Should print size with units (B, KB, MB, or GB)
  expect_message(tdr_cache_info(mock_tdr), "Total size:")
})

test_that("tdr_cache_info skips NULL slot manifests", {
  cache_dir <- file.path(tempdir(), "tdr_test_info_null_slot")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  
  m1 <- .tdr_cache_write("data", cache_dir, "clustering.ids", "s1")
  
  mock_tdr <- list(
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
  mock_tdr <- list(
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