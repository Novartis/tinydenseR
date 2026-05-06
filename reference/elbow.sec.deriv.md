# Find elbow point using second derivative method (internal)

Detects the "elbow" in a curve by finding where the second derivative
(rate of change of the slope) is maximized. This identifies the point
where the curve transitions from steep to flat, commonly used for
determining optimal dimensionality in PCA or Laplacian Eigenmaps.

## Usage

``` r
elbow.sec.deriv(x, smooth = TRUE, df = NULL, sort.order = c("asc", "desc"))
```

## Arguments

- x:

  Numeric vector of values (e.g., eigenvalues, explained variance). Must
  have at least 4 values for second derivative calculation.

- smooth:

  Logical, whether to apply smoothing spline before computing second
  derivative. Default TRUE. Only applied if length(x) \> 6.

- df:

  Degrees of freedom for smoothing spline. If NULL (default),
  auto-calculated as min(length(x) \* 0.7, 10). Lower values = smoother
  fit.

- sort.order:

  Character, either "asc" or "desc". Specifies the expected ordering of
  meaningful values:

  - `"asc"` - values increase with importance (Laplacian Eigenmaps)

  - `"desc"` - values decrease with importance (PCA eigenvalues)

  The algorithm works on sorted data, so this ensures correct
  interpretation.

## Value

List with three elements:

- `index`:

  Integer index of elbow point in original `x` vector

- `value`:

  Numeric value of `x` at the elbow point

- `sec.deriv`:

  Maximum second derivative value (diagnostic)

## Details

The algorithm:

1.  Sorts values according to `sort.order`

2.  Optionally smooths with spline to reduce noise

3.  Normalizes to `[0,1]` range for scale-invariance

4.  Computes second derivative using `diff(..., differences = 2)`

5.  Returns index where second derivative is largest

Smoothing is recommended for noisy data but requires at least 7 points.
The `df` parameter controls smoothing strength: lower values = smoother.
