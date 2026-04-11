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
#
# Cache files are written to the system temporary directory and are
# intended to be ephemeral — they exist only for the duration of the
# current R session.  An explicit cleanup finalizer ensures the cache
# is removed when the session ends (see .tdr_cache_register_cleanup).
# ──────────────────────────────────────────────────────────────────────

# Current schema version — bump when cache format changes
.TDR_CACHE_SCHEMA_VERSION <- 3L

# ── Package-private environment for session-end cleanup ──────────────
# A strong-referenced environment that lives for the duration of the R
# session.  reg.finalizer(..., onexit = TRUE) is attached to it once
# per session so that all registered cache directories are removed when
# the session exits (or when the environment is garbage-collected).
.tdr_cleanup_env <- new.env(parent = emptyenv())
.tdr_cleanup_env$registered_dirs <- character(0)
.tdr_cleanup_env$finalizer_set   <- FALSE

#' Register a cache directory for automatic removal at session end
#'
#' Adds \code{cache_dir} to a package-private registry and, on first
#' call, attaches a \code{reg.finalizer(..., onexit = TRUE)} to a
#' package-private environment so that all registered directories are
#' deleted when the R session terminates.
#'
#' @param cache_dir Character path of the cache directory to register.
#' @return Invisible \code{NULL}.
#' @keywords internal
.tdr_cache_register_cleanup <- function(cache_dir) {
  .tdr_cleanup_env$registered_dirs <- unique(c(
    .tdr_cleanup_env$registered_dirs, cache_dir
  ))
  if (!isTRUE(.tdr_cleanup_env$finalizer_set)) {
    reg.finalizer(
      e = .tdr_cleanup_env,
      f = function(env) {
        dirs <- env$registered_dirs
        for (d in dirs) {
          if (dir.exists(d)) {
            unlink(d, recursive = TRUE, force = TRUE)
          }
        }
      },
      onexit = TRUE
    )
    .tdr_cleanup_env$finalizer_set <- TRUE
  }
  invisible(NULL)
}

#' Remove a single cache directory from the cleanup registry
#'
#' Called by \code{tdr_cache_cleanup()} after manually deleting a cache
#' so that the finalizer does not attempt a redundant \code{unlink()}.
#'
#' @param cache_dir Character path to deregister.
#' @return Invisible \code{NULL}.
#' @keywords internal
.tdr_cache_deregister <- function(cache_dir) {
  .tdr_cleanup_env$registered_dirs <- setdiff(
    .tdr_cleanup_env$registered_dirs, cache_dir
  )
  invisible(NULL)
}

#' Return the root temporary cache directory for this R session
#'
#' All \code{get.map()} cache runs are stored under
#' \code{<tempdir()>/tinydenseR_cache/}.  The directory is created
#' lazily on first call.
#'
#' @return Character scalar — the session-level cache root path.
#' @keywords internal
.tdr_cache_root <- function() {
  root <- file.path(tempdir(), "tinydenseR_cache")
  if (!dir.exists(root)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  root
}

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
#' @param slot_name Character – one of \code{"clustering"},
#'   \code{"celltyping"}, \code{"nearest.lm"},
#'   \code{"fuzzy.graphs"}.
#' @param sample_name Character – sample identifier (used as file stem).
#' @param compress Passed to \code{saveRDS}. Default \code{FALSE} for speed.
#' @return An attributed path string: a character scalar with
#'   \code{attr(, "schema_v")} and \code{attr(, "bytes")} set.
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
  
  # Return attributed path string
  path <- target
  attr(path, "schema_v") <- .TDR_CACHE_SCHEMA_VERSION
  attr(path, "bytes")    <- file.size(target)
  path
}

#' Read a cached slot from disk
#'
#' Accepts either a new-style attributed path string (schema v3+) or a
#' legacy metadata list (schema v2) with a \code{$path} element.
#'
#' @param meta_record An attributed path string (character scalar with
#'   \code{attr(, "schema_v")} and \code{attr(, "bytes")}), OR a
#'   legacy metadata list with \code{$path} and \code{$schema_v}.
#' @return The deserialized R object.
#' @keywords internal
.tdr_cache_read <- function(meta_record) {
  # Detect format: path string vs legacy metadata list
  if (is.character(meta_record) && length(meta_record) == 1L) {
    path     <- meta_record
    schema_v <- attr(meta_record, "schema_v")
  } else if (is.list(meta_record) && !is.null(meta_record$path)) {
    # Legacy v2 metadata list
    path     <- meta_record$path
    schema_v <- meta_record$schema_v
  } else {
    stop("Invalid cache record: expected an attributed path string or a ",
         "metadata list with $path.")
  }
  
  if (!file.exists(path)) {
    stop("Cache file missing: ", path,
         "\nRe-run get.map() to regenerate or set .cache.on.disk = FALSE.")
  }
  
  # Schema version gate
  if (!is.null(schema_v) &&
      schema_v != .TDR_CACHE_SCHEMA_VERSION) {
    stop("Cache schema mismatch: file has v", schema_v,
         " but current code expects v", .TDR_CACHE_SCHEMA_VERSION,
         ".\nRe-run get.map() to regenerate.")
  }
  
  readRDS(file = path)
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
#' per-cell data from \code{@cellmap}.  Supports the unified cellmap
#' structure with four top-level slots: \code{clustering},
#' \code{celltyping}, \code{nearest.lm}, and \code{fuzzy.graphs}.
#' Clustering and celltyping have an inner \code{$ids} sub-list;
#' nearest.lm and fuzzy.graphs are flat named lists.
#'
#' Auto-detects on-disk path strings (\code{is.character(val) &&
#' length(val) == 1 && file.exists(val)}) and loads them with
#' \code{readRDS}.
#'
#' @param .tdr.obj tinydenseR object.
#' @param slot_name Character – one of \code{"clustering"},
#'   \code{"celltyping"}, \code{"nearest.lm"}, \code{"fuzzy.graphs"}.
#'   Legacy names (\code{"clustering.ids"}, \code{"celltyping.ids"},
#'   \code{"nearest.landmarks"}, \code{"fuzzy.graph"}) are mapped
#'   automatically.
#' @param sample_name Character – sample identifier.
#' @return The R object for that sample/slot.
#' @keywords internal
.tdr_get_map_slot <- function(.tdr.obj, slot_name, sample_name) {
  
  # Legacy name mapping
  slot_name <- switch(slot_name,
    "clustering.ids"    = "clustering",
    "celltyping.ids"    = "celltyping",
    "nearest.landmarks" = "nearest.lm",
    "fuzzy.graph"       = "fuzzy.graphs",
    slot_name
  )
  
  # Navigate to the correct sub-list
  val <- switch(slot_name,
    "clustering"   = .tdr.obj@cellmap$clustering$ids[[sample_name]],
    "celltyping"   = .tdr.obj@cellmap$celltyping$ids[[sample_name]],
    "nearest.lm"   = .tdr.obj@cellmap$nearest.lm[[sample_name]],
    "fuzzy.graphs" = .tdr.obj@cellmap$fuzzy.graphs[[sample_name]],
    stop("Unknown slot_name: ", slot_name)
  )
  
  if (is.null(val)) {
    stop("No entry for slot '", slot_name,
         "', sample '", sample_name, "'.")
  }
  
  # Auto-detect on-disk path string and lazy-load
  if (is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))) {
    if (!file.exists(val)) {
      stop("Cache file missing: ", val)
    }
    return(.tdr_cache_read(val))
  }
  if (is.character(val) && length(val) == 1L && file.exists(val)) {
    return(.tdr_cache_read(val))
  }
  
  val
}

#' Access a full slot as a named list (all samples), reading from cache lazily
#'
#' Returns a named list over all samples for the given slot, loading each
#' element from disk when the entry is an on-disk path string.
#'
#' @param .tdr.obj tinydenseR object.
#' @param slot_name Character – one of \code{"clustering"},
#'   \code{"celltyping"}, \code{"nearest.lm"}, \code{"fuzzy.graphs"}.
#'   Legacy names are mapped automatically.
#' @return A named list of R objects (one per sample).
#' @keywords internal
.tdr_get_map_slot_all <- function(.tdr.obj, slot_name) {
  
  # Legacy name mapping
  slot_name <- switch(slot_name,
    "clustering.ids"    = "clustering",
    "celltyping.ids"    = "celltyping",
    "nearest.landmarks" = "nearest.lm",
    "fuzzy.graph"       = "fuzzy.graphs",
    slot_name
  )
  
  entries <- switch(slot_name,
    "clustering"   = .tdr.obj@cellmap$clustering$ids,
    "celltyping"   = .tdr.obj@cellmap$celltyping$ids,
    "nearest.lm"   = .tdr.obj@cellmap$nearest.lm,
    "fuzzy.graphs" = .tdr.obj@cellmap$fuzzy.graphs,
    stop("Unknown slot_name: ", slot_name)
  )
  
  if (is.null(entries)) return(NULL)
  
  # Lazy-load any on-disk path strings
  lapply(entries, function(val) {
    if (is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))) {
      if (!file.exists(val)) stop("Cache file missing: ", val)
      .tdr_cache_read(val)
    } else if (is.character(val) && length(val) == 1L && file.exists(val)) {
      .tdr_cache_read(val)
    } else {
      val
    }
  })
}

#' Access per-cell map data for a single sample
#'
#' Retrieves per-cell data (cluster IDs, cell type IDs, nearest landmarks,
#' or fuzzy graph) for a single sample, transparently reading from on-disk
#' cache when caching is active or from in-memory slots otherwise.
#'
#' @param x A \code{\linkS4class{TDRObj}} processed through \code{get.map()}.
#' @param .slot Character. One of \code{"clustering"},
#'   \code{"celltyping"}, \code{"nearest.lm"}, \code{"fuzzy.graphs"}.
#'   Legacy names (\code{"clustering.ids"}, \code{"celltyping.ids"},
#'   \code{"nearest.landmarks"}, \code{"fuzzy.graph"}) are accepted for
#'   backward compatibility.
#' @param .sample Character. Sample identifier (must match a name in
#'   \code{names(x@@cells)}).
#' @return The per-cell object for that sample and slot: a named character
#'   vector (for \code{*ids}), an integer matrix (for
#'   \code{nearest.lm}), or a sparse matrix (for
#'   \code{fuzzy.graphs}).
#'
#' @seealso \code{\link{get.map}}, \code{\link{tdr_cache_validate}}
#'
#' @examples
#' \dontrun{
#' # Get cell type labels for the first sample
#' ct <- get.cellmap(x = lm.obj,
#'                   .slot = "celltyping",
#'                   .sample = names(lm.obj$cells)[1])
#' }
#' @export
get.cellmap <- function(x, .slot, .sample) {
  if (!is.TDRObj(x)) stop("x must be a TDRObj.")
  .tdr_get_map_slot(x, slot_name = .slot, sample_name = .sample)
}

#' Validate that all on-disk cache files are intact
#'
#' Walks the \code{@cellmap} sub-slots for path strings and checks file
#' existence.  This function is called automatically (in quiet mode) when
#' entering \code{get.lm()}, \code{get.pbDE()}, and \code{get.plsD()} so
#' that broken caches are caught early.
#'
#' Because the cache is ephemeral (stored under \code{tempdir()} and
#' cleaned up at session end), heavyweight checksum verification is no
#' longer performed.  Only file existence and schema version are checked.
#'
#' @param .tdr.obj A tinydenseR object.
#' @param .verbose Logical; if \code{TRUE}, print a summary.  Default
#'   \code{TRUE}.
#' @return The (unmodified) \code{.tdr.obj}, invisibly.  Stops with an
#'   error if any cached file is missing.
#' @export
tdr_cache_validate <- function(.tdr.obj,
                               .verbose = TRUE) {
  
  # Collect all path strings from @cellmap sub-slots
  .is_path_string <- function(val) {
    is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))
  }
  
  # Define which sub-structures to walk
  slot_lists <- list(
    clustering   = .tdr.obj@cellmap$clustering$ids,
    celltyping   = .tdr.obj@cellmap$celltyping$ids,
    nearest.lm   = .tdr.obj@cellmap$nearest.lm,
    fuzzy.graphs = .tdr.obj@cellmap$fuzzy.graphs
  )
  
  n_ok      <- 0L
  n_missing <- 0L
  problems  <- character(0)
  has_paths <- FALSE
  
  for (slot_name in names(slot_lists)) {
    entries <- slot_lists[[slot_name]]
    if (is.null(entries)) next
    for (sample_name in names(entries)) {
      val <- entries[[sample_name]]
      if (is.null(val)) next
      if (!.is_path_string(val)) next
      
      has_paths <- TRUE
      
      # Schema version check
      sv <- attr(val, "schema_v")
      if (!is.null(sv) && sv != .TDR_CACHE_SCHEMA_VERSION) {
        stop("Cache schema mismatch: file has v", sv,
             " but current code expects v", .TDR_CACHE_SCHEMA_VERSION,
             ".\nRe-run get.map() to regenerate.")
      }
      
      if (!file.exists(val)) {
        n_missing <- n_missing + 1L
        problems <- c(problems,
                      paste0("MISSING: ", slot_name, "/", sample_name,
                             " -> ", val))
      } else {
        n_ok <- n_ok + 1L
      }
    }
  }
  
  if (!has_paths) {
    if (isTRUE(.verbose)) message("No on-disk cache active; nothing to validate.")
    return(invisible(.tdr.obj))
  }
  
  if (length(problems) > 0L) {
    stop("Cache validation failed (\n  ",
         paste(problems, collapse = "\n  "),
         "\n)\nRe-run get.map() to regenerate, or run ",
         "tdr_cache_cleanup() and re-map.")
  }
  
  if (isTRUE(.verbose)) {
    msg <- paste0("Cache OK: ", n_ok, " file(s) verified.")
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
  # Quick check: any path strings at all?
  has_cache <- !is.null(.tdr.obj@config$.cache.root) ||
    any(vapply(
      list(.tdr.obj@cellmap$clustering$ids,
           .tdr.obj@cellmap$celltyping$ids,
           .tdr.obj@cellmap$nearest.lm,
           .tdr.obj@cellmap$fuzzy.graphs),
      function(lst) {
        if (is.null(lst)) return(FALSE)
        any(vapply(lst, function(v) {
          is.character(v) && length(v) == 1L && !is.null(attr(v, "schema_v"))
        }, logical(1)))
      },
      logical(1)
    ))
  if (!has_cache) return(invisible(NULL))
  tdr_cache_validate(.tdr.obj, .verbose = FALSE)
  invisible(NULL)
}

#' Remove all cached files for a tinydenseR object
#'
#' Deletes the on-disk cache directory, deregisters it from the
#' session-end cleanup finalizer, and clears the path strings from
#' \code{@cellmap}.  After cleanup, the cellmap entries will be
#' \code{NULL} and must be regenerated via \code{get.map()}.
#'
#' @param .tdr.obj A tinydenseR object.
#' @return Updated \code{.tdr.obj} with cache removed.
#' @export
tdr_cache_cleanup <- function(.tdr.obj) {
  # Derive cache root from @config$.cache.root or from path strings
  cache_root <- .tdr.obj@config$.cache.root
  if (is.null(cache_root)) {
    # Fallback: derive from any path string (dirname/dirname → run dir)
    for (lst in list(.tdr.obj@cellmap$clustering$ids,
                     .tdr.obj@cellmap$nearest.lm,
                     .tdr.obj@cellmap$fuzzy.graphs,
                     .tdr.obj@cellmap$celltyping$ids)) {
      if (is.null(lst)) next
      for (val in lst) {
        if (is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))) {
          cache_root <- dirname(dirname(val))
          break
        }
      }
      if (!is.null(cache_root)) break
    }
  }
  
  if (!is.null(cache_root) && dir.exists(cache_root)) {
    unlink(cache_root, recursive = TRUE)
    .tdr_cache_deregister(cache_root)
    message("Removed cache directory: ", cache_root)
  }
  
  # Clear path strings from cellmap
  .tdr.obj@cellmap$clustering$ids  <- NULL
  .tdr.obj@cellmap$celltyping$ids  <- NULL
  .tdr.obj@cellmap$nearest.lm      <- NULL
  .tdr.obj@cellmap$fuzzy.graphs    <- NULL
  .tdr.obj@config$.cache.root      <- NULL
  .tdr.obj
}


#' Print a human-readable summary of the on-disk cache state
#'
#' Reports whether on-disk caching is active, the cache directory, schema
#' version, number of cached slots and files, and total size on disk.
#' Useful for interactive inspection and debugging.
#'
#' @param .tdr.obj A tinydenseR object.
#' @return A list (invisible) with components \code{active}, \code{root},
#'   \code{schema_v}, \code{slots} (character vector of cached slot names),
#'   \code{n_samples}, \code{n_files}, and \code{total_bytes}.
#' @examples
#' \dontrun{
#' # After running get.map() with on-disk caching
#' tdr_cache_info(lm.cells)
#' # On-disk cache: ACTIVE
#' # Directory:     /tmp/RtmpXXXXXX/tinydenseR_cache/run_20260302_143012_a7f3c1b2
#' # Schema:        v3
#' # Slots:         clustering, celltyping, nearest.lm, fuzzy.graphs
#' # Samples:       6
#' # Files:         24  (4 slots x 6 samples)
#' # Total size:    142.3 MB
#' }
#' @seealso \code{\link{tdr_cache_validate}}, \code{\link{tdr_cache_cleanup}},
#'   \code{\link{get.map}}
#' @export
tdr_cache_info <- function(.tdr.obj) {
  
  .is_path_string <- function(val) {
    is.character(val) && length(val) == 1L && !is.null(attr(val, "schema_v"))
  }
  
  slot_lists <- list(
    clustering   = .tdr.obj@cellmap$clustering$ids,
    celltyping   = .tdr.obj@cellmap$celltyping$ids,
    nearest.lm   = .tdr.obj@cellmap$nearest.lm,
    fuzzy.graphs = .tdr.obj@cellmap$fuzzy.graphs
  )
  
  # Count files and accumulate bytes
  n_files     <- 0L
  total_bytes <- 0
  sample_set  <- character(0)
  active_slots <- character(0)
  schema_v    <- NULL
  
  for (sn in names(slot_lists)) {
    entries <- slot_lists[[sn]]
    if (is.null(entries)) next
    slot_has_paths <- FALSE
    for (sm in names(entries)) {
      val <- entries[[sm]]
      if (is.null(val) || !.is_path_string(val)) next
      slot_has_paths <- TRUE
      n_files <- n_files + 1L
      sample_set <- union(sample_set, sm)
      b <- attr(val, "bytes")
      if (!is.null(b)) total_bytes <- total_bytes + b
      if (is.null(schema_v)) schema_v <- attr(val, "schema_v")
    }
    if (slot_has_paths) active_slots <- c(active_slots, sn)
  }
  
  if (n_files == 0L) {
    message("On-disk cache: INACTIVE")
    message("(Re-run get.map(.cache.on.disk = TRUE) to enable.)")
    return(invisible(list(active = FALSE)))
  }
  
  n_samples <- length(sample_set)
  cache_root <- .tdr.obj@config$.cache.root
  
  # Human-readable size
  fmt_size <- function(b) {
    if (b < 1024)         return(paste0(b, " B"))
    if (b < 1024^2)       return(sprintf("%.1f KB", b / 1024))
    if (b < 1024^3)       return(sprintf("%.1f MB", b / 1024^2))
    sprintf("%.1f GB", b / 1024^3)
  }
  
  message("On-disk cache: ACTIVE")
  if (!is.null(cache_root)) message("Directory:     ", cache_root)
  if (!is.null(schema_v))   message("Schema:        v", schema_v)
  message("Slots:         ", paste(active_slots, collapse = ", "))
  message("Samples:       ", n_samples)
  message("Files:         ", n_files,
          "  (", length(active_slots), " slot(s) x ", n_samples, " sample(s))")
  message("Total size:    ", fmt_size(total_bytes))
  
  invisible(list(
    active      = TRUE,
    root        = cache_root,
    schema_v    = schema_v,
    slots       = active_slots,
    n_samples   = n_samples,
    n_files     = n_files,
    total_bytes = total_bytes
  ))
}
