#' Find elbow point using second derivative method (internal)
#'
#' Detects the "elbow" in a curve by finding where the second derivative
#' (rate of change of the slope) is maximized. This identifies the point
#' where the curve transitions from steep to flat, commonly used for
#' determining optimal dimensionality in PCA or Laplacian Eigenmaps.
#' 
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Sorts values according to \code{sort.order}
#'   \item Optionally smooths with spline to reduce noise
#'   \item Normalizes to \verb{[0,1]} range for scale-invariance
#'   \item Computes second derivative using \code{diff(..., differences = 2)}
#'   \item Returns index where second derivative is largest
#' }
#' 
#' Smoothing is recommended for noisy data but requires at least 7 points.
#' The \code{df} parameter controls smoothing strength: lower values = smoother.
#' 
#' @param x Numeric vector of values (e.g., eigenvalues, explained variance).
#'   Must have at least 4 values for second derivative calculation.
#' @param smooth Logical, whether to apply smoothing spline before computing
#'   second derivative. Default TRUE. Only applied if length(x) > 6.
#' @param df Degrees of freedom for smoothing spline. If NULL (default),
#'   auto-calculated as min(length(x) * 0.7, 10). Lower values = smoother fit.
#' @param sort.order Character, either "asc" or "desc". Specifies the expected
#'   ordering of meaningful values:
#'   \itemize{
#'     \item \code{"asc"} - values increase with importance (Laplacian Eigenmaps)
#'     \item \code{"desc"} - values decrease with importance (PCA eigenvalues)
#'   }
#'   The algorithm works on sorted data, so this ensures correct interpretation.
#' 
#' @return List with three elements:
#'   \describe{
#'     \item{\code{index}}{Integer index of elbow point in original \code{x} vector}
#'     \item{\code{value}}{Numeric value of \code{x} at the elbow point}
#'     \item{\code{sec.deriv}}{Maximum second derivative value (diagnostic)}
#'   }
#' 
#' @keywords internal
#' 
elbow.sec.deriv <- 
  function(x,
           smooth = TRUE, 
           df = NULL, 
           sort.order = c("asc", "desc")) {
    
    sort.order <- 
      match.arg(arg = sort.order)
    
    if(length(x = x) < 4) {
      stop("Cannot compute elbow point: need at least 4 values for second derivative calculation.\n",
           "Current length: ", length(x))
    }
    
    orig.x <- x
    ord <- 
      order(x, 
            decreasing = (sort.order == "desc"))
    x_sorted <- x[ord]
    
    # Apply smoothing to reduce noise (only if enough data points)
    if(smooth &&
       length(x = x_sorted) > 6) {
      
      # Auto-calculate df: 70% of length, capped at 10 for stability
      if(is.null(x = df)) df <- min(length(x = x_sorted) * 0.7, 10)
      
      # Clamp df to valid range: [2, n-1] required by smooth.spline
      df <- max(2,
                min(df, 
                    length(x = x_sorted) - 1))
      
      x_sorted <- 
        stats::smooth.spline(x = seq_along(along.with = x_sorted),
                             y = x_sorted, 
                             df = df)$y
    }
    
    
    # Normalize to [0,1] for scale-invariant comparison
    rng <- max(x_sorted) - min(x_sorted)
    if (rng == 0) {
      stop("Cannot detect elbow: all values are identical (flat sequence).\n",
           "This indicates no variation in the data.")
    }
    
    norm.vals <-
      (x_sorted - min(x_sorted)) / 
      rng
    
    # Compute second derivative: measures curvature (rate of slope change)
    sec.deriv <- 
      diff(x = norm.vals, 
           differences = 2)
    
    # Find maximum curvature point and map back to original indices
    elbow.idx.sorted <- which.max(x = sec.deriv) + 1  # +1 accounts for diff() reducing length
    elbow.idx.orig <- ord[elbow.idx.sorted]
    
    list(
      index = elbow.idx.orig,
      value = orig.x[elbow.idx.orig],
      sec.deriv = sec.deriv[which.max(x = sec.deriv)]
    )
  }