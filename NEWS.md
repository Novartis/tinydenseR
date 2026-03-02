# tinydenseR 0.0.2.0001

## On-disk caching for `get.map()` (new feature)

`get.map()` now defaults to **on-disk caching** (`.cache.on.disk = TRUE`),
serializing four large per-sample slots to disk instead of holding them
in memory:

- `clustering$ids`
- `celltyping$ids`
- `nearest.landmarks`
- `fuzzy.graph`

This dramatically reduces peak RAM usage for large single-cell datasets.
All downstream functions (`goi.summary()`, `get.lm()`, `get.pbDE()`,
`get.specDE()`, `get.nmfDE()`, `get.plsDE()`) read from the cache
transparently.

### Cache management functions

| Function | Purpose |
|---|---|
| `tdr_cache_info()` | Print a human-readable summary of cache state (directory, size, slots, samples) |
| `tdr_cache_validate()` | Verify that all cached files exist (and optionally check MD5 checksums with `.verify.checksum = TRUE`) |
| `tdr_cache_cleanup()` | Delete the cache directory and strip metadata from the object |

### Hardening

- **MD5 checksums** are stored per file at write time and can be verified
  with `tdr_cache_validate(.verify.checksum = TRUE)` or
  `.tdr_cache_read(meta, verify = TRUE)`.
- **Schema versioning** — cache metadata carries a `schema_v` field; reads
  error immediately when a version mismatch is detected.
- **Crash-safe writes** — files are written to a `.rds.tmp` temp file and
  atomically renamed; orphaned temp files are cleaned up automatically.
- **File-level locking** — when the `filelock` package is installed
  (added to `Suggests`), a lock is acquired during the atomic rename to
  prevent concurrent-write clobbering.
- **Automatic validation** — `get.lm()`, `get.pbDE()`, and
  `.validate.DE.inputs()` (used by `get.specDE()`, `get.nmfDE()`,
  `get.plsDE()`) silently validate the cache on entry.

### Disabling caching

Pass `.cache.on.disk = FALSE` to `get.map()` to keep all data in memory
(original behaviour).
