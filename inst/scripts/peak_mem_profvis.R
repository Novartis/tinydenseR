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

# Compute peak memory and elapsed time from a profvis profile
peak_mem_profvis <- function(p,
                             input  = c("auto", "bytes", "MiB", "MB"),
                             output = c("MiB", "MB", "bytes"),
                             verbose = FALSE,
                             peak_mode = c("absolute", "over_baseline")) {
  input     <- match.arg(input)
  output    <- match.arg(output)
  peak_mode <- match.arg(peak_mode)
  
  # utility
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  # -------------------------------------------------------------------------
  # 1) Extract prof table
  # -------------------------------------------------------------------------
  prof <- NULL
  is_widget <- inherits(p, "htmlwidget")
  if (is_widget) {
    prof <- p$x$message$prof %||% p$x$prof
    if (is.null(prof)) stop("Could not locate 'prof' table inside profvis object.")
  } else if (is.data.frame(p)) {
    prof <- p
  } else {
    stop("Argument 'p' must be a profvis htmlwidget or a data.frame.")
  }
  
  # -------------------------------------------------------------------------
  # 2) Per-sample memory deltas
  # -------------------------------------------------------------------------
  if (!("meminc" %in% names(prof)) && !("memalloc" %in% names(prof))) {
    stop("No memory columns found. Expected 'meminc' or 'memalloc'.")
  }
  deltas_raw <- if ("meminc" %in% names(prof)) prof$meminc else c(0, diff(prof$memalloc))
  deltas_raw <- deltas_raw[is.finite(deltas_raw)]
  if (!length(deltas_raw)) stop("Memory delta series is empty after filtering.")
  
  # -------------------------------------------------------------------------
  # 3) Auto-detect input memory unit (bytes vs MiB)
  # -------------------------------------------------------------------------
  unit_guess  <- NULL
  diagnostics <- list()
  if (input == "auto") {
    nz        <- deltas_raw[abs(deltas_raw) > 0]
    max_abs   <- if (length(nz)) max(abs(nz)) else 0
    frac_prop <- if (length(nz)) mean(abs(nz - round(nz)) > 1e-6) else 0
    # Try peaks under both interpretations
    peak_if_bytes_bytes <- diff(range(cumsum(deltas_raw)))          # treat as BYTES
    peak_if_mib_bytes   <- diff(range(cumsum(deltas_raw) * 1024^2)) # treat as MiB->bytes
    peak_if_bytes_mib   <- peak_if_bytes_bytes / 1024^2
    peak_if_mib_mib     <- peak_if_mib_bytes / 1024^2
    if (frac_prop > 0.10 && max_abs < 1024) {
      unit_guess <- "MiB"
    } else if (max_abs >= 8 * 1024^2) {
      unit_guess <- "bytes"
    } else {
      plausible_bytes <- peak_if_bytes_mib > 0.1 && peak_if_bytes_mib < 64 * 1024^2
      plausible_mib   <- peak_if_mib_mib   > 0.1 && peak_if_mib_mib   < 64 * 1024^2
      if (plausible_mib && !plausible_bytes) {
        unit_guess <- "MiB"
      } else if (plausible_bytes && !plausible_mib) {
        unit_guess <- "bytes"
      } else {
        unit_guess <- "MiB"
        diagnostics$ambiguous <- TRUE
      }
    }
    diagnostics$max_abs_raw       <- max_abs
    diagnostics$frac_non_integer  <- frac_prop
    diagnostics$peak_if_bytes_mib <- peak_if_bytes_mib
    diagnostics$peak_if_mib_mib   <- peak_if_mib_mib
    diagnostics$unit_guess        <- unit_guess
  } else {
    unit_guess <- input
  }
  if (unit_guess %in% c("MiB", "MB") && input == "bytes") {
    warning("Input appears pre-scaled (MiB/MB), but 'input=\"bytes\"' was supplied.")
  }
  
  # Convert to BYTES
  deltas_bytes <- switch(unit_guess,
                         bytes = deltas_raw,
                         MiB   = deltas_raw * 1024^2,
                         MB    = deltas_raw * 1e6)
  cum_bytes <- cumsum(deltas_bytes)
  
  # Peak memory
  peak_bytes <- if (peak_mode == "absolute") max(cum_bytes, na.rm = TRUE)
  else                          diff(range(cum_bytes))
  peak_out   <- switch(output,
                       MiB   = peak_bytes / 1024^2,
                       MB    = peak_bytes / 1e6,
                       bytes = peak_bytes)
  
  # -------------------------------------------------------------------------
  # 4) Elapsed time from sampling interval — robust extraction & normalization
  # -------------------------------------------------------------------------
  interval_candidates_raw <- c()
  timing_mode <- NULL
  if (is_widget) {
    interval_candidates_raw <- c(interval_candidates_raw,
                                 p$x$message$interval,
                                 p$x$options$interval,
                                 p$x$interval)
    timing_mode <- p$x$options$timing %||% p$x$message$timing %||% "elapsed"
  }
  interval_candidates_raw <- as.numeric(interval_candidates_raw)
  interval_candidates_raw <- interval_candidates_raw[is.finite(interval_candidates_raw) & interval_candidates_raw > 0]
  
  normalize_interval <- function(v) {
    # profvis UI commonly stores 10 (ms) while message uses 0.01 (s)
    if (v >= 1 && v <= 1000 && abs(v - round(v)) < 1e-9) {
      # looks like an integer number of milliseconds (e.g., 10)
      return(v / 1000)
    }
    v
  }
  interval_candidates_s <- vapply(interval_candidates_raw, normalize_interval, numeric(1))
  
  # Heuristic: prefer the smallest normalized positive value (typical sampling ~0.005-0.05s)
  interval_s <- if (length(interval_candidates_s)) min(interval_candidates_s) else NA_real_
  if (!is.finite(interval_s)) interval_s <- 0.01
  
  # Count top-level samples; if depth is missing (e.g., synthetic df), this returns 0
  n_samples_depth1 <- if ("depth" %in% names(prof)) sum(prof$depth == 1L, na.rm = TRUE) else 0L
  elapsed_s  <- n_samples_depth1 * interval_s
  elapsed_min <- elapsed_s / 60
  
  timing <- list(
    interval_candidates_raw = interval_candidates_raw,
    interval_candidates_s   = interval_candidates_s,
    interval_s              = interval_s,
    n_samples_depth1        = n_samples_depth1,
    elapsed_s               = elapsed_s,
    elapsed_min             = elapsed_min,
    timing_mode             = timing_mode
  )
  
  diagnostics$input_param    <- input
  diagnostics$output_param   <- output
  diagnostics$unit_used      <- unit_guess
  diagnostics$samples        <- nrow(prof)
  diagnostics$nonzero_ratio  <- mean(deltas_raw != 0)
  diagnostics$peak_out       <- peak_out
  diagnostics$peak_mode      <- peak_mode
  
  if (verbose) {
    cat("\n[peak_mem_profvis] Diagnostics\n",
        "--------------------------------\n",
        sprintf("Samples (rows in prof):  %d\n", diagnostics$samples),
        sprintf("Top-level samples:       %d\n", n_samples_depth1),
        sprintf("Timing mode:             %s\n", timing_mode %||% "unknown"),
        sprintf("Interval candidates raw: %s\n",
                if (length(interval_candidates_raw)) paste(interval_candidates_raw, collapse = ", ") else "none"),
        sprintf("Interval candidates (s): %s\n",
                if (length(interval_candidates_s)) paste(sprintf("%.4f", interval_candidates_s), collapse = ", ") else "none"),
        sprintf("Chosen interval (s):     %.4f\n", interval_s),
        sprintf("Elapsed:                 %.3f s (%.3f min)\n",
                elapsed_s, elapsed_min),
        sprintf("Non-zero delta ratio:    %.2f\n", diagnostics$nonzero_ratio),
        sprintf("Raw max |delta|:         %.3f (raw units)\n", diagnostics$max_abs_raw %||% NA_real_),
        sprintf("Frac non-integer deltas: %.2f\n", diagnostics$frac_non_integer %||% NA_real_),
        sprintf("Peak if BYTES:           %.3f MiB\n", diagnostics$peak_if_bytes_mib %||% NA_real_),
        sprintf("Peak if MiB:             %.3f MiB\n", diagnostics$peak_if_mib_mib %||% NA_real_),
        sprintf("Unit guess:              %s\n", diagnostics$unit_guess),
        sprintf("Peak mode:               %s\n", peak_mode),
        sprintf("Output unit:             %s\n", output),
        sprintf("Final peak:              %.3f %s\n\n", peak_out, output),
        sep = "")
    if (!is.null(diagnostics$ambiguous) && diagnostics$ambiguous) {
      warning("Unit detection was ambiguous; defaulting to MiB. ",
              "Consider specifying input = \"bytes\" or input = \"MiB\" explicitly.")
    }
  }
  
  structure(list(
    peak             = peak_out,
    unit             = output,
    peak_bytes       = peak_bytes,
    cum_series_bytes = cum_bytes,
    timing           = timing,
    diagnostics      = diagnostics
  ), class = "peak_mem_profvis")
}

# Print method (unchanged)
print.peak_mem_profvis <- function(x, ...) {
  cat(sprintf("Peak memory: %.3f %s | Elapsed: %.3f s (%.3f min)\n",
              x$peak, x$unit, x$timing$elapsed_s, x$timing$elapsed_min))
  invisible(x)
}