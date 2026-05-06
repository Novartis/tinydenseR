# Access fuzzy density matrices from a TDRObj

Returns one of the three fuzzy density layers stored in `@density` after
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
has been run.

## Usage

``` r
get.density(x, ...)

# S3 method for class 'TDRObj'
get.density(x, .which = c("raw", "norm", "log.norm", "fdens", "Y"), ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  (or Seurat / SCE / HDF5AnnData wrapping one via
  [`GetTDR`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md)).

- ...:

  Additional arguments passed to methods.

- .which:

  Character(1). Which density layer to return:

  `"raw"`

  :   Pre-size-factor-normalization fuzzy density sums (landmarks ×
      samples). Each entry is \\\sum_c F_j\[c,l\]\\, the total fuzzy
      edge weight connecting sample \\j\\'s cells to landmark \\l\\.

  `"norm"`

  :   Size-factor-normalized fuzzy density (landmarks × samples). Equals
      `raw / size.factors` (column-wise). This is the primary analytical
      layer. Aliases: `"fdens"` (deprecated).

  `"log.norm"`

  :   Log2-transformed normalized density: `log2(norm + 0.5)`. Used by
      [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
      and
      [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md).
      Aliases: `"Y"` (deprecated).

## Value

A numeric matrix (landmarks × samples), or `NULL` if the requested layer
is not available.

## Details

The three density layers form a deterministic chain computed by
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md):

\$\$\mathrm{raw}\[l,j\] = \sum\_{c \in j} F_j\[c,l\]\$\$

where \\F_j\\ is the fuzzy graph (UMAP edge weights) for sample \\j\\,
\\c\\ indexes cells, and \\l\\ indexes landmarks.

\$\$\mathrm{norm} = t(t(\mathrm{raw}) \\/\\ \mathrm{size.factors})\$\$

with \\\mathrm{size.factors}\[j\] = n_j / \bar{n}\\ (cells per sample
divided by the mean cell count across samples).

\$\$\mathrm{log.norm} = \log_2(\mathrm{norm} + 0.5)\$\$

The pseudo-count of 0.5 stabilizes landmarks with zero density.

## Deprecated aliases

`"fdens"` maps to `"norm"` and `"Y"` maps to `"log.norm"`. These aliases
are retained for backward compatibility with code written against
tinydenseR \< 0.0.3.

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
(producer),
[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`get.embedding`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md),
[`as.SummarizedExperiment.TDRObj`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.TDRObj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running the pipeline:
raw_dens  <- get.density(tdr, "raw")       # pre-normalization
norm_dens <- get.density(tdr, "norm")       # size-factor normalized
log_dens  <- get.density(tdr, "log.norm")   # log2(norm + 0.5)

# Verify invariants:
all.equal(norm_dens, t(t(raw_dens) / tdr@density$size.factors))
all.equal(log_dens, log2(norm_dens + 0.5))
} # }
```
