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

library(testthat)
library(tinydenseR)

# ── Shared fixture: minimal cyto pipeline through get.map() ──
# Uses bimodal data (two well-separated blobs) to guarantee 2+ clusters

.build_mapped_obj <- function(seed = 42, cache_on_disk = FALSE, n_samples = 4) {
  
  set.seed(seed)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  blob2 <- matrix(rnorm(n * m, mean = 5, sd = 0.5), nrow = n, ncol = m)
  
  .cells <- lapply(seq_len(n_samples), function(i) {
    mat <- rbind(blob1 + rnorm(1, 0, 0.1), blob2 + rnorm(1, 0, 0.1))
    dimnames(mat) <- list(
      paste0("s", i, "_", seq_len(nrow(mat))),
      paste0("M", seq_len(m))
    )
    uri <- tempfile(fileext = ".RDS")
    saveRDS(object = mat, file = uri, compress = FALSE)
    return(uri)
  })
  names(.cells) <- paste0("sample", seq_len(n_samples))
  
  .meta <- data.frame(
    row.names = names(.cells),
    group = rep(c("A", "B"), length.out = n_samples)
  )
  
  obj <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("M", seq_len(m)),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(
      .k = 5,
      .scale = FALSE,
      .verbose = FALSE,
      .seed = seed
    ) |>
    get.map(
      .verbose = FALSE,
      .seed = seed,
      .cache.on.disk = cache_on_disk
    )
  
  obj
}

# Build a celltyping map from whatever clusters the object has
.make_ct_map <- function(obj) {
  cls <- levels(obj@landmark.annot$clustering$ids)
  n <- length(cls)
  if (n < 2) {
    # If only 1 cluster, map it to a single celltype
    return(list("TypeA" = cls))
  }
  half <- floor(n / 2)
  list(
    "TypeA" = cls[seq_len(half)],
    "TypeB" = cls[(half + 1):n]
  )
}

# ──────────────────────────────────────────────────────────────────────
# U1: celltyping() runs after get.map() without error
# ──────────────────────────────────────────────────────────────────────

test_that("U1: celltyping() runs after get.map() without error", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  expect_no_error(
    obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  )
  expect_true(is.TDRObj(obj2))
})

# ──────────────────────────────────────────────────────────────────────
# U2: landmark.annot$celltyping$ids is correctly updated after late celltyping
# ──────────────────────────────────────────────────────────────────────

test_that("U2: celltyping IDs are correctly updated after late celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  ids <- obj2@landmark.annot$celltyping$ids
  expect_s3_class(ids, "factor")
  expect_true(all(ids %in% c("TypeA", "TypeB")))
  expect_equal(length(ids), nrow(obj2@assay$expr))
})

# ──────────────────────────────────────────────────────────────────────
# U3: median.exprs and pheatmap are recomputed
# ──────────────────────────────────────────────────────────────────────

test_that("U3: median.exprs and pheatmap are recomputed after late celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  expect_true(!is.null(obj2@results$celltyping$median.exprs))
  expect_true(is.matrix(obj2@results$celltyping$median.exprs))
  expect_equal(sort(rownames(obj2@results$celltyping$median.exprs)),
               sort(c("TypeA", "TypeB")))
  
  expect_true(!is.null(obj2@results$celltyping$pheatmap))
  expect_s3_class(obj2@results$celltyping$pheatmap, "pheatmap")
})

# ──────────────────────────────────────────────────────────────────────
# U4: celltyping() before get.map() still works (backward compat)
# ──────────────────────────────────────────────────────────────────────

test_that("U4: celltyping() before get.map() still works", {
  set.seed(42)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, 0, 0.5), n, m)
  blob2 <- matrix(rnorm(n * m, 5, 0.5), n, m)
  mat1 <- rbind(blob1, blob2)
  dimnames(mat1) <- list(paste0("s1_", 1:nrow(mat1)), paste0("M", 1:m))
  f1 <- tempfile(fileext = ".RDS")
  saveRDS(mat1, f1, compress = FALSE)
  .cells <- list(sample1 = f1)
  .meta <- data.frame(row.names = "sample1", group = "A")
  
  obj <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("M", 1:m),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 42) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = 42)
  
  ct_map <- .make_ct_map(obj)
  
  expect_no_error(
    obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  )
  
  # fdens should still be NULL since get.map was not run
  expect_null(obj2@density$fdens)
  
  # But landmark-level celltyping should be set
  expect_true(!is.null(obj2@landmark.annot$celltyping$ids))
  expect_true(!is.null(obj2@results$celltyping$median.exprs))
})

# ──────────────────────────────────────────────────────────────────────
# U5: .refresh_celltyping() skips get.map-dependent slots when fdens is NULL
# ──────────────────────────────────────────────────────────────────────

test_that("U5: .refresh_celltyping skips map-dependent slots when fdens is NULL", {
  set.seed(42)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, 0, 0.5), n, m)
  blob2 <- matrix(rnorm(n * m, 5, 0.5), n, m)
  mat1 <- rbind(blob1, blob2)
  dimnames(mat1) <- list(paste0("s1_", 1:nrow(mat1)), paste0("M", 1:m))
  f1 <- tempfile(fileext = ".RDS")
  saveRDS(mat1, f1, compress = FALSE)
  .cells <- list(sample1 = f1)
  .meta <- data.frame(row.names = "sample1", group = "A")
  
  obj <- setup.tdr.obj(
    .cells = .cells,
    .meta = .meta,
    .markers = paste0("M", 1:m),
    .assay.type = "cyto",
    .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = 42) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = 42)
  
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  # cellmap$celltype.ids should be NULL (no get.map run)
  expect_null(obj2@cellmap$celltype.ids)
  
  # composition celltyping should not exist
  expect_null(obj2@density$composition$celltyping)
})

# ──────────────────────────────────────────────────────────────────────
# I1: Full late-celltyping pipeline — verify cellmap$celltype.ids
#     reflects updated labels
# ──────────────────────────────────────────────────────────────────────

test_that("I1: cellmap$celltype.ids reflects updated labels after late celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  # Cell-level IDs should now be populated with the new celltype labels
  ct_ids <- .tdr_get_map_slot_all(obj2, "celltyping.ids")
  expect_true(!is.null(ct_ids))
  
  all_labels <- unique(unlist(lapply(ct_ids, unique)))
  # Should contain TypeA, TypeB, and possibly ..low.confidence..
  expect_true(any(c("TypeA", "TypeB") %in% all_labels))
})

# ──────────────────────────────────────────────────────────────────────
# I1b: Composition matrices are refreshed
# ──────────────────────────────────────────────────────────────────────

test_that("I1b: composition matrices are refreshed after late celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  cell_count <- obj2@density$composition$celltyping$cell.count
  cell_perc  <- obj2@density$composition$celltyping$cell.perc
  
  expect_true(!is.null(cell_count))
  expect_true(!is.null(cell_perc))
  expect_true(is.matrix(cell_count))
  expect_true(is.matrix(cell_perc))
  
  # Column names should be the celltype labels (possibly including ..low.confidence..)
  expect_true(any(c("TypeA", "TypeB") %in% colnames(cell_count)))
  
  # Percentages should sum to ~100
  expect_true(all(abs(rowSums(cell_perc) - 100) < 1e-6))
})

# ──────────────────────────────────────────────────────────────────────
# N1: fdens and Y are unchanged after celltyping re-run
# ──────────────────────────────────────────────────────────────────────

test_that("N1: fdens and Y are not changed by late celltyping", {
  obj <- .build_mapped_obj()
  
  fdens_before <- obj@density$fdens
  Y_before     <- obj@density$Y
  
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  expect_identical(obj2@density$fdens, fdens_before)
  expect_identical(obj2@density$Y, Y_before)
})

# ──────────────────────────────────────────────────────────────────────
# N2: clustering, graphs, assay, landmark.embed are unchanged
# ──────────────────────────────────────────────────────────────────────

test_that("N2: non-celltyping slots are unchanged after late celltyping", {
  obj <- .build_mapped_obj()
  
  clustering_before  <- obj@landmark.annot$clustering
  graphs_before      <- obj@graphs
  assay_before       <- obj@assay
  embed_before       <- obj@landmark.embed
  cluster_ids_before <- .tdr_get_map_slot_all(obj, "clustering.ids")
  
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  expect_identical(obj2@landmark.annot$clustering, clustering_before)
  expect_identical(obj2@graphs, graphs_before)
  expect_identical(obj2@assay, assay_before)
  expect_identical(obj2@landmark.embed, embed_before)
  
  cluster_ids_after <- .tdr_get_map_slot_all(obj2, "clustering.ids")
  expect_identical(cluster_ids_after, cluster_ids_before)
})

# ──────────────────────────────────────────────────────────────────────
# .celltyping.map is stored for provenance
# ──────────────────────────────────────────────────────────────────────

test_that("celltyping.map and mode are stored for provenance", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  expect_identical(obj2@landmark.annot$celltyping$map, ct_map)
  expect_identical(obj2@landmark.annot$celltyping$mode, "cluster_map")
})

# ──────────────────────────────────────────────────────────────────────
# On-disk cache path: celltyping IDs are rewritten
# ──────────────────────────────────────────────────────────────────────

test_that("I4: on-disk cached celltyping IDs are rewritten after late celltyping", {
  obj <- .build_mapped_obj(cache_on_disk = TRUE)
  
  cache <- obj@density$.cache
  skip_if(is.null(cache) || !isTRUE(cache$on.disk),
          "On-disk caching not active")
  
  ct_map <- .make_ct_map(obj)
  obj2 <- celltyping(obj, ct_map, .verbose = FALSE)
  
  # Verify cache manifests are updated
  cache2 <- obj2@density$.cache
  expect_true(!is.null(cache2$manifests$celltyping.ids))
  
  # Read back from disk and verify labels
  for (sn in names(obj2@cells)) {
    ids <- .tdr_get_map_slot(obj2, "celltyping.ids", sn)
    expect_true(any(c("TypeA", "TypeB") %in% unique(ids)))
  }
})

# ──────────────────────────────────────────────────────────────────────
# I2: get.lm() trad$celltyping is invalidated with warning after late celltyping
# ──────────────────────────────────────────────────────────────────────

test_that("I2: trad$celltyping is invalidated after late celltyping", {
  obj <- .build_mapped_obj()
  ct_map <- .make_ct_map(obj)
  
  # First celltyping + get.lm
  obj <- celltyping(obj, ct_map, .verbose = FALSE)
  
  .design <- model.matrix(~ group, data = obj@metadata)
  
  suppressWarnings(
    obj <- get.lm(obj, .design = .design, .verbose = FALSE)
  )
  
  # Verify trad$celltyping fit exists
  expect_true(!is.null(obj@results$lm[["default"]]$trad$celltyping$fit))
  
  # Re-run celltyping → should invalidate trad$celltyping with a warning
  expect_warning(
    obj2 <- celltyping(obj, ct_map, .verbose = FALSE),
    "Invalidated trad\\$celltyping"
  )
  
  expect_null(obj2@results$lm[["default"]]$trad$celltyping)
})

# ──────────────────────────────────────────────────────────────────────
# R2: Results match between old (celltyping → get.map) and new 
#     (get.map → celltyping) workflows
# ──────────────────────────────────────────────────────────────────────

test_that("R2: old vs new workflow produces equivalent composition", {
  seed <- 42
  set.seed(seed)
  n <- 200; m <- 5
  blob1 <- matrix(rnorm(n * m, 0, 0.5), n, m)
  blob2 <- matrix(rnorm(n * m, 5, 0.5), n, m)
  
  # Common data — 2 samples
  mats <- list(
    sample1 = {
      mat <- rbind(blob1, blob2)
      dimnames(mat) <- list(paste0("s1_", 1:nrow(mat)), paste0("M", 1:m))
      mat
    },
    sample2 = {
      mat <- rbind(blob1 + 0.1, blob2 + 0.1)
      dimnames(mat) <- list(paste0("s2_", 1:nrow(mat)), paste0("M", 1:m))
      mat
    }
  )
  
  save_cells <- function(mats) {
    lapply(mats, function(x) {
      uri <- tempfile(fileext = ".RDS")
      saveRDS(object = x, file = uri, compress = FALSE)
      return(uri)
    })
  }
  
  .meta <- data.frame(
    row.names = c("sample1", "sample2"),
    group = c("A", "B")
  )
  
  # Old workflow: celltyping() → get.map()
  .cells1 <- save_cells(mats)
  obj_old <- setup.tdr.obj(
    .cells = .cells1, .meta = .meta,
    .markers = paste0("M", 1:m),
    .assay.type = "cyto", .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = seed)
  
  ct_map <- .make_ct_map(obj_old)
  
  obj_old <- obj_old |>
    celltyping(ct_map, .verbose = FALSE) |>
    get.map(.verbose = FALSE, .seed = seed)
  
  # New workflow: get.map() → celltyping()
  .cells2 <- save_cells(mats)
  obj_new <- setup.tdr.obj(
    .cells = .cells2, .meta = .meta,
    .markers = paste0("M", 1:m),
    .assay.type = "cyto", .verbose = FALSE
  ) |>
    get.landmarks(.verbose = FALSE, .seed = seed) |>
    get.graph(.k = 5, .scale = FALSE, .verbose = FALSE, .seed = seed) |>
    get.map(.verbose = FALSE, .seed = seed) |>
    celltyping(ct_map, .verbose = FALSE)
  
  # Compare cell percentages (should be identical since the fuzzy graph is
  # the same and the vote logic is the same)
  old_perc <- obj_old@density$composition$celltyping$cell.perc
  new_perc <- obj_new@density$composition$celltyping$cell.perc
  
  # Align columns (order may differ)
  common_cols <- intersect(colnames(old_perc), colnames(new_perc))
  expect_true(length(common_cols) > 0)
  
  expect_equal(
    old_perc[, common_cols, drop = FALSE],
    new_perc[, common_cols, drop = FALSE],
    tolerance = 1e-10
  )
})

# ──────────────────────────────────────────────────────────────────────
# label.confidence is stored in config by get.map()
# ──────────────────────────────────────────────────────────────────────

test_that("label.confidence is stored in config by get.map", {
  obj <- .build_mapped_obj()
  expect_equal(obj@config$label.confidence, 0.8)
})
