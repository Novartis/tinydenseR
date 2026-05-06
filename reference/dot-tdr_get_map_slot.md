# Access a per-sample map slot, transparently reading from cache if on-disk

This is the single entry point downstream code should use to access
per-cell data from `@cellmap`. Supports the unified cellmap structure
with four top-level slots: `clustering`, `celltyping`, `nearest.lm`, and
`fuzzy.graphs`. Clustering and celltyping have an inner `$ids` sub-list;
nearest.lm and fuzzy.graphs are flat named lists.

## Usage

``` r
.tdr_get_map_slot(.tdr.obj, slot_name, sample_name)
```

## Arguments

- .tdr.obj:

  tinydenseR object.

- slot_name:

  Character – one of `"clustering"`, `"celltyping"`, `"nearest.lm"`,
  `"fuzzy.graphs"`. Legacy names (`"clustering.ids"`,
  `"celltyping.ids"`, `"nearest.landmarks"`, `"fuzzy.graph"`) are mapped
  automatically.

- sample_name:

  Character – sample identifier.

## Value

The R object for that sample/slot.

## Details

Auto-detects on-disk path strings
(`is.character(val) && length(val) == 1 && file.exists(val)`) and loads
them with `readRDS`.
