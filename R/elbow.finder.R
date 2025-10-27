#' Find elbow point using second derivative method
#' 
#' @param x Numeric vector (e.g., PCA eigenvalues in decreasing order)
#' @param smooth Logical, whether to apply smoothing spline
#' @param df Degrees of freedom for smoothing (auto-calculated if NULL)
#' @param sort.order Character, either `asc` or `desc` for sorting order
#'   - `asc`  for Laplacian Eigenmaps eigenvalues (increasing).
#'   - `desc` for PCA eigenvalues (decreasing).
#' @return List with index, value, and second derivative at elbow point
#' 
elbow.sec.deriv <- 
  function(x,
           smooth = TRUE, 
           df = NULL, 
           sort.order = c("asc", "desc")) {
    
    sort.order <- 
      match.arg(arg = sort.order)
    
    if(length(x = x) < 4) stop("Need at least 4 values for second derivative")
    
    orig.x <- x
    ord <- 
      order(x, 
            decreasing = (sort.order == "desc"))
    x_sorted <- x[ord]
    
    if(smooth &&
       length(x = x_sorted) > 6) {
      
      if(is.null(x = df)) df <- min(length(x = x_sorted) * 0.7, 10)
      
      # clamp df into valid range for smooth.spline
      df <- max(2,
                min(df, 
                    length(x = x_sorted) - 1))
      
      x_sorted <- 
        stats::smooth.spline(x = seq_along(along.with = x_sorted),
                             y = x_sorted, 
                             df = df)$y
    }
    
    
    rng <- max(x_sorted) - min(x_sorted)
    if (rng == 0) stop("Sequence is flat; no elbow detectable")
    
    norm.vals <-
      (x_sorted - min(x_sorted)) / 
      rng
    
    sec.deriv <- 
      diff(x = norm.vals, 
           differences = 2)
    
    elbow.idx.sorted <- which.max(x = sec.deriv) + 1
    elbow.idx.orig <- ord[elbow.idx.sorted]
    
    list(
      index = elbow.idx.orig,
      value = orig.x[elbow.idx.orig],
      sec.deriv = sec.deriv[which.max(x = sec.deriv)]
    )
  }