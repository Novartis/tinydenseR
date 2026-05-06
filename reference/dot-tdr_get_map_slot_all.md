# Access a full slot as a named list (all samples), reading from cache lazily

Returns a named list over all samples for the given slot, loading each
element from disk when the entry is an on-disk path string.

## Usage

``` r
.tdr_get_map_slot_all(.tdr.obj, slot_name)
```

## Arguments

- .tdr.obj:

  tinydenseR object.

- slot_name:

  Character – one of `"clustering"`, `"celltyping"`, `"nearest.lm"`,
  `"fuzzy.graphs"`. Legacy names are mapped automatically.

## Value

A named list of R objects (one per sample).
