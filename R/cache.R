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
# Internal on-disk caching helpers for large get.map() slots
# ──────────────────────────────────────────────────────────────────────

# Current schema version — bump when cache format changes
.TDR_CACHE_SCHEMA_VERSION <- 1L

#' Generate a unique run key for cache directory isolation
#'
#' Creates a key from timestamp + a random hex string to avoid collisions
#' across concurrent runs and across datasets.
#'
#' @return Character scalar, e.g. "run_20260302_143012_a7f3c1b2"
#' @keywords internal
.tdr_make_run_key <- function() {
  paste0("run_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
         paste0(sprintf("%02x", sample.int(256, 4) - 1L), collapse = ""))
}

#' Write one cached slot to disk (crash-safe via atomic rename)
#'
#' Writes \code{object} to an RDS file inside \code{cache_dir/slot_name/}
#' using an atomic rename pattern (write to a \code{.rds.tmp} file, then
#' \code{file.rename()} to the final \code{.rds}).  When the \pkg{filelock}
#' package is installed, a shared file lock is held during the rename to
#' prevent concurrent writes from clobbering each other.
#'
#' @param object R object to serialize.
#' @param cache_dir Root cache directory for the current run.
#' @param slot_name Character – one of \code{"clustering.ids"},
#'   \code{"celltyping.ids"}, \code{"nearest.landmarks"},
#'   \code{"fuzzy.graph"}.
#' @param sample_name Character – sample identifier (used as file stem).
#' @param compress Passed to \code{saveRDS}. Default \code{FALSE} for speed.
#' @return A metadata list with path, slot, sample, bytes, md5, and schema_v.
#' @keywords internal
.tdr_cache_write <- function(object,
                             cache_dir,
                             slot_name,
                             sample_name,
                             compress = FALSE) {
  
  slot_dir <- file.path(cache_dir, slot_name)
  dir.create(slot_dir, recursive = TRUE, showWarnings = FALSE)
  
  target <- file.path(slot_dir, paste0(sample_name, ".rds"))
  
  # Atomic write: save to temp file then rename
  tmp <- tempfile(tmpdir = slot_dir, fileext = ".rds.tmp")
  saveRDS(object = object, file = tmp, compress = compress)
  
  # Optional file-level locking for concurrent-write safety
  lck <- NULL
  if (requireNamespace("filelock", quietly = TRUE)) {
    lck <- filelock::lock(path = paste0(target, ".lock"), timeout = 10000)
  }
  on.exit({
    if (!is.null(lck)) filelock::unlock(lck)
    # Remove lock file (best-effort)
    lockfile <- paste0(target, ".lock")
    if (file.exists(lockfile)) try(file.remove(lockfile), silent = TRUE)
  }, add = TRUE)
  
  file.rename(from = tmp, to = target)
  
  # Compute md5 checksum for later verification
  md5 <- unname(tools::md5sum(target))
  
  list(
    path     = target,
    slot     = slot_name,
    sample   = sample_name,
    bytes    = file.size(target),
    md5      = md5,
    schema_v = .TDR_CACHE_SCHEMA_VERSION
  )
}

#' Read a cached slot from disk
#'
#' @param meta_record A metadata list produced by \code{.tdr_cache_write}.
#' @param verify Logical; if \code{TRUE}, verify the MD5 checksum before
#'   deserializing.  Default \code{FALSE} for speed.
#' @return The deserialized R object.
#' @keywords internal
.tdr_cache_read <- function(meta_record, verify = FALSE) {
  if (!file.exists(meta_record$path)) {
    stop("Cache file missing: ", meta_record$path,
         "\nRe-run get.map() to regenerate or set .cache.on.disk = FALSE.")
  }
  
  # Schema version gate
  if (!is.null(meta_record$schema_v) &&
      meta_record$schema_v != .TDR_CACHE_SCHEMA_VERSION) {
    stop("Cache schema mismatch: file has v", meta_record$schema_v,
         " but current code expects v", .TDR_CACHE_SCHEMA_VERSION,
         ".\nRe-run get.map() to regenerate.")
  }
  
  # Optional checksum verification
  if (isTRUE(verify) && !is.null(meta_record$md5)) {
    current_md5 <- unname(tools::md5sum(meta_record$path))
    if (!identical(current_md5, meta_record$md5)) {
      stop("Cache checksum mismatch for ", meta_record$path,
           "\nExpected: ", meta_record$md5,
           "\nGot:      ", current_md5,
           "\nThe file may be corrupted. Re-run get.map().")
    }
  }
  
  readRDS(file = meta_record$path)
}

#' Remove orphaned .rds.tmp files left by interrupted writes
#'
#' @param cache_dir Cache root directory to sweep.
#' @keywords internal
.tdr_cache_sweep_orphans <- function(cache_dir) {
  if (!dir.exists(cache_dir)) return(invisible(NULL))
  tmps <- list.files(cache_dir, pattern = "\\.rds\\.tmp$",
                     recursive = TRUE, full.names = TRUE)
  if (length(tmps) > 0L) {
    file.remove(tmps)
    message("Cleaned up ", length(tmps), " orphaned temp file(s).")
  }
  invisible(NULL)
}

#' Access a per-sample map slot, transparently reading from cache if on-disk
#'
#' This is the single entry point downstream code should use to access
#' \code{$map$clustering$ids[[sample]]}, \code{$map$celltyping$ids[[sample]]},
#' \code{$map$nearest.landmarks[[sample]]}, and
#' \code{$map$fuzzy.graph[[sample]]}.
#'
#' @param .tdr.obj tinydenseR object.
#' @param slot_name Character – one of "clustering.ids", "celltyping.ids",
#'   "nearest.landmarks", "fuzzy.graph".
#' @param sample_name Character – sample identifier.
#' @return The R object for that sample/slot.
#' @keywords internal
.tdr_get_map_slot <- function(.tdr.obj, slot_name, sample_name) {
  
  cache <- .tdr.obj@map$.cache
  
  if (!is.null(cache) && isTRUE(cache$on.disk)) {
    # Read from disk
    meta <- cache$manifests[[slot_name]][[sample_name]]
    if (is.null(meta)) {
      stop("No cache entry for slot '", slot_name,
           "', sample '", sample_name, "'.")
    }
    return(.tdr_cache_read(meta))
  }
  
  # In-memory fallback
  switch(slot_name,
         "clustering.ids"    = .tdr.obj@map$clustering$ids[[sample_name]],
         "celltyping.ids"    = .tdr.obj@map$celltyping$ids[[sample_name]],
         "nearest.landmarks" = .tdr.obj@map$nearest.landmarks[[sample_name]],
         "fuzzy.graph"       = .tdr.obj@map$fuzzy.graph[[sample_name]],
         stop("Unknown slot_name: ", slot_name))
}

#' Access a full slot as a named list (all samples), reading from cache lazily
#'
#' Returns a named list over all samples for the given slot, loading each
#' element from disk when caching is active.
#'
#' @param .tdr.obj tinydenseR object.
#' @param slot_name Character – one of "clustering.ids", "celltyping.ids",
#'   "nearest.landmarks", "fuzzy.graph".
#' @return A named list of R objects (one per sample).
#' @keywords internal
.tdr_get_map_slot_all <- function(.tdr.obj, slot_name) {
  
  cache <- .tdr.obj@map$.cache
  
  if (!is.null(cache) && isTRUE(cache$on.disk)) {
    manifests <- cache$manifests[[slot_name]]
    if (is.null(manifests)) return(NULL)
    lapply(X = manifests, FUN = .tdr_cache_read)
  } else {
    switch(slot_name,
           "clustering.ids"    = .tdr.obj@map$clustering$ids,
           "celltyping.ids"    = .tdr.obj@map$celltyping$ids,
           "nearest.landmarks" = .tdr.obj@map$nearest.landmarks,
           "fuzzy.graph"       = .tdr.obj@map$fuzzy.graph,
           stop("Unknown slot_name: ", slot_name))
  }
}

#' Validate that all on-disk cache files are intact
#'
#' Checks every entry in the cache manifest for file existence and,
#' optionally, MD5 checksum integrity.  This function is called
#' automatically (in quiet mode) when entering \code{get.lm()},
#' \code{get.pbDE()}, \code{get.specDE()}, \code{get.nmfDE()}, and
#' \code{get.plsDE()} so that broken caches are caught early.
#'
#' @param .tdr.obj A tinydenseR object with \code{$map$.cache} populated.
#' @param .verify.checksum Logical; if \code{TRUE}, also verify MD5
#'   checksums.  Default \code{FALSE} for speed.
#' @param .verbose Logical; if \code{TRUE}, print a summary.  Default
#'   \code{TRUE}.
#' @return The (unmodified) \code{.tdr.obj}, invisibly.  Stops with an
#'   error if any file is missing or any checksum mismatches and
#'   \code{.verify.checksum = TRUE}.
#' @export
tdr_cache_validate <- function(.tdr.obj,
                               .verify.checksum = FALSE,
                               .verbose = TRUE) {
  
  cache <- .tdr.obj@map$.cache
  
  if (is.null(cache) || !isTRUE(cache$on.disk)) {
    if (isTRUE(.verbose)) message("No on-disk cache active; nothing to validate.")
    return(invisible(.tdr.obj))
  }
  
  # Schema version check
  if (!is.null(cache$schema_v) &&
      cache$schema_v != .TDR_CACHE_SCHEMA_VERSION) {
    stop("Cache schema mismatch: manifest has v", cache$schema_v,
         " but current code expects v", .TDR_CACHE_SCHEMA_VERSION,
         ".\nRe-run get.map() to regenerate.")
  }
  
  n_ok      <- 0L
  n_missing <- 0L
  n_bad_md5 <- 0L
  problems  <- character(0)
  
  for (slot_name in names(cache$manifests)) {
    for (sample_name in names(cache$manifests[[slot_name]])) {
      meta <- cache$manifests[[slot_name]][[sample_name]]
      if (is.null(meta)) next
      
      if (!file.exists(meta$path)) {
        n_missing <- n_missing + 1L
        problems <- c(problems,
                      paste0("MISSING: ", slot_name, "/", sample_name,
                             " -> ", meta$path))
      } else {
        n_ok <- n_ok + 1L
        if (isTRUE(.verify.checksum) && !is.null(meta$md5)) {
          current_md5 <- unname(tools::md5sum(meta$path))
          if (!identical(current_md5, meta$md5)) {
            n_bad_md5 <- n_bad_md5 + 1L
            problems <- c(problems,
                          paste0("CHECKSUM MISMATCH: ", slot_name, "/",
                                 sample_name, " (expected ", meta$md5,
                                 ", got ", current_md5, ")"))
          }
        }
      }
    }
  }
  
  if (length(problems) > 0L) {
    stop("Cache validation failed (\n  ",
         paste(problems, collapse = "\n  "),
         "\n)\nRe-run get.map() to regenerate, or run ",
         "tdr_cache_cleanup() and re-map.")
  }
  
  if (isTRUE(.verbose)) {
    msg <- paste0("Cache OK: ", n_ok, " file(s) verified")
    if (isTRUE(.verify.checksum)) msg <- paste0(msg, " (with MD5)")
    msg <- paste0(msg, ".")
    message(msg)
  }
  
  invisible(.tdr.obj)
}

#' Quiet internal cache validation for DE entry points
#'
#' Runs the same checks as \code{tdr_cache_validate()} but silently;
#' it only stops on errors and never emits messages.  Designed to be
#' injected at the top of \code{get.lm()}, \code{get.pbDE()}, etc.
#'
#' @param .tdr.obj A tinydenseR object.
#' @return \code{NULL} invisibly.
#' @keywords internal
.tdr_cache_validate_quiet <- function(.tdr.obj) {
  cache <- .tdr.obj@map$.cache
  if (is.null(cache) || !isTRUE(cache$on.disk)) return(invisible(NULL))
  tdr_cache_validate(.tdr.obj, .verify.checksum = FALSE, .verbose = FALSE)
  invisible(NULL)
}

#' Remove all cached files for a tinydenseR object
#'
#' Deletes the on-disk cache directory and strips cache metadata from the
#' object. After cleanup, the large slots (clustering$ids, celltyping$ids,
#' nearest.landmarks, fuzzy.graph) will be \code{NULL} and must be regenerated
#' via \code{get.map()}.
#'
#' @param .tdr.obj A tinydenseR object with \code{$map$.cache} populated.
#' @return Updated \code{.tdr.obj} with cache removed.
#' @export
tdr_cache_cleanup <- function(.tdr.obj) {
  cache_root <- .tdr.obj@map$.cache$root
  if (!is.null(cache_root) && dir.exists(cache_root)) {
    unlink(cache_root, recursive = TRUE)
    message("Removed cache directory: ", cache_root)
  }
  .tdr.obj@map$.cache <- NULL
  .tdr.obj
}


#' Print a human-readable summary of the on-disk cache state
#'
#' Reports whether on-disk caching is active, the cache directory, schema
#' version, number of cached slots and files, and total size on disk.
#' Useful for interactive inspection and debugging.
#'
#' @param .tdr.obj A tinydenseR object (with or without \code{$map$.cache}).
#' @return A list (invisible) with components \code{active}, \code{root},
#'   \code{schema_v}, \code{slots} (character vector of cached slot names),
#'   \code{n_samples}, \code{n_files}, and \code{total_bytes}.
#' @examples
#' \dontrun{
#' # After running get.map() with on-disk caching
#' tdr_cache_info(lm.cells)
#' # On-disk cache: ACTIVE
#' # Directory:     /data/tdr_cache/run_20260302_143012_a7f3c1b2
#' # Schema:        v1
#' # Slots:         clustering.ids, celltyping.ids, nearest.landmarks, fuzzy.graph
#' # Samples:       6
#' # Files:         24  (4 slots x 6 samples)
#' # Total size:    142.3 MB
#' }
#' @seealso \code{\link{tdr_cache_validate}}, \code{\link{tdr_cache_cleanup}},
#'   \code{\link{get.map}}
#' @export
tdr_cache_info <- function(.tdr.obj) {
  
  cache <- .tdr.obj@map$.cache
  
  if (is.null(cache) || !isTRUE(cache$on.disk)) {
    message("On-disk cache: INACTIVE")
    message("(Re-run get.map(.cache.on.disk = TRUE) to enable.)")
    return(invisible(list(active = FALSE)))
  }
  
  slot_names <- names(cache$manifests)
  slot_names <- slot_names[!vapply(cache$manifests[slot_names], is.null, logical(1))]
  
  # Count files and accumulate bytes
  n_files    <- 0L
  total_bytes <- 0
  sample_set <- character(0)
  
  for (sn in slot_names) {
    entries <- cache$manifests[[sn]]
    for (sm in names(entries)) {
      meta <- entries[[sm]]
      if (is.null(meta)) next
      n_files <- n_files + 1L
      sample_set <- union(sample_set, sm)
      if (!is.null(meta$bytes)) {
        total_bytes <- total_bytes + meta$bytes
      }
    }
  }
  
  n_samples <- length(sample_set)
  
  # Human-readable size
  fmt_size <- function(b) {
    if (b < 1024)         return(paste0(b, " B"))
    if (b < 1024^2)       return(sprintf("%.1f KB", b / 1024))
    if (b < 1024^3)       return(sprintf("%.1f MB", b / 1024^2))
    sprintf("%.1f GB", b / 1024^3)
  }
  
  message("On-disk cache: ACTIVE")
  message("Directory:     ", cache$root)
  message("Schema:        v", cache$schema_v)
  message("Slots:         ", paste(slot_names, collapse = ", "))
  message("Samples:       ", n_samples)
  message("Files:         ", n_files,
          "  (", length(slot_names), " slot(s) x ", n_samples, " sample(s))")
  message("Total size:    ", fmt_size(total_bytes))
  
  invisible(list(
    active      = TRUE,
    root        = cache$root,
    schema_v    = cache$schema_v,
    slots       = slot_names,
    n_samples   = n_samples,
    n_files     = n_files,
    total_bytes = total_bytes
  ))
}
