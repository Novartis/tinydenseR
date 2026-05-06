# Convert a TDRObj to SummarizedExperiment

Converts a tinydenseR
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
into a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Usage

``` r
# S3 method for class 'TDRObj'
as.SummarizedExperiment(x, ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  (must have
  [`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  completed).

- ...:

  Additional arguments (currently unused).

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Details

Rows represent tinydenseR **landmarks** (not genes/proteins). The assays
stored are:

- `counts`:

  Raw fuzzy graph density sums before size-factor normalization, from
  `@density$raw`.

- `normcounts`:

  Size-factor-normalized fuzzy density (landmarks × samples), from
  `@density$norm`.

- `logcounts`:

  log2(normcounts + 0.5), used by
  [`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  for linear modeling, from `@density$log.norm`.

`rowData` contains all stored clustering and celltyping solutions.
`colData` contains sample-level metadata. The full TDRObj is preserved
in `metadata(se)$tdr.obj` for downstream tinydenseR analysis via
[`GetTDR`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md).

## See also

[`GetTDR`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md),
[`RunTDR`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert TDRObj to SummarizedExperiment
se <- as.SummarizedExperiment(tdr.obj)

# Access density data
SummarizedExperiment::assay(se, "logcounts")[1:5, 1:5]

# Round-trip: extract TDRObj for full tinydenseR analysis
tdr <- GetTDR(se)
} # }
```
