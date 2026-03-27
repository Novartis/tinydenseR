library(testthat)
library(tinydenseR)

.make_seurat5_multilayer_fixture <- function(split_sample = FALSE) {
  counts <- Matrix::Matrix(
    data = matrix(
      c(
        1, 0, 3, 0, 5, 0,
        0, 2, 0, 4, 0, 6,
        7, 8, 9, 10, 11, 12
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE,
    dimnames = list(
      paste0("g", 1:3),
      paste0("c", 1:6)
    )
  )

  sample_id <- c("s1", "s1", "s2", "s2", "s3", "s3")

  meta <- data.frame(
    sample_id = sample_id,
    row.names = colnames(counts)
  )

  seu <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  seu[["RNA"]] <- methods::as(seu[["RNA"]], "Assay5")

  if (split_sample) {
    SeuratObject::LayerData(seu, assay = "RNA", layer = "counts.a") <- counts[, c("c1", "c3", "c5")]
    SeuratObject::LayerData(seu, assay = "RNA", layer = "counts.b") <- counts[, c("c2", "c4", "c6")]
    seu@meta.data$layer_id <- c("a", "b", "a", "b", "a", "b")
  } else {
    SeuratObject::LayerData(seu, assay = "RNA", layer = "counts.s12") <- counts[, c("c1", "c2", "c3", "c4")]
    SeuratObject::LayerData(seu, assay = "RNA", layer = "counts.s3") <- counts[, c("c5", "c6")]
    seu@meta.data$layer_id <- c("s12", "s12", "s12", "s12", "s3", "s3")
  }

  list(seu = seu, counts = counts)
}

test_that("get.cells.Seurat5 returns one file per sample and respects .meta row order", {
  skip_if_not_installed("SeuratObject")

  fixture <- .make_seurat5_multilayer_fixture(split_sample = FALSE)

  sample_meta <- data.frame(
    group = c("A", "B", "C"),
    row.names = c("s3", "s1", "s2")
  )

  cells <- tinydenseR::get.cells.Seurat5(
    .seurat.obj = fixture$seu,
    .meta = sample_meta,
    .sample.var = "sample_id",
    .assay = "RNA",
    .layer.pattern = "counts.",
    .layer.name.source = "layer_id",
    .min.cells.per.sample = 2,
    .verbose = FALSE
  )

  on.exit(unlink(unlist(cells)), add = TRUE)

  expect_identical(names(cells), rownames(sample_meta))
  expect_equal(length(cells), 3)

  expect_true(all(file.exists(unlist(cells))))

  mats <- lapply(cells, readRDS)
  expect_true(all(vapply(mats, methods::is, logical(1), "dgCMatrix")))

  expect_equal(dim(mats$s1), c(3, 2))
  expect_equal(dim(mats$s2), c(3, 2))
  expect_equal(dim(mats$s3), c(3, 2))

  expect_setequal(colnames(mats$s1), c("c1", "c2"))
  expect_setequal(colnames(mats$s2), c("c3", "c4"))
  expect_setequal(colnames(mats$s3), c("c5", "c6"))
})

test_that("get.cells.Seurat5 errors when a sample has cells in multiple layers", {
  skip_if_not_installed("SeuratObject")

  fixture <- .make_seurat5_multilayer_fixture(split_sample = TRUE)

  sample_meta <- data.frame(
    group = c("A", "B", "C"),
    row.names = c("s1", "s2", "s3")
  )

  expect_error(
    suppressWarnings(
      tinydenseR::get.cells.Seurat5(
        .seurat.obj = fixture$seu,
        .meta = sample_meta,
        .sample.var = "sample_id",
        .assay = "RNA",
        .layer.pattern = "counts.",
        .layer.name.source = NULL,
        .min.cells.per.sample = 2,
        .verbose = FALSE
      )
    ),
    "across multiple layers"
  )
})

test_that("get.cells.Seurat5 validates .layer.name.source column", {
  skip_if_not_installed("SeuratObject")

  fixture <- .make_seurat5_multilayer_fixture(split_sample = FALSE)

  sample_meta <- data.frame(
    group = c("A", "B", "C"),
    row.names = c("s1", "s2", "s3")
  )

  expect_error(
    tinydenseR::get.cells.Seurat5(
      .seurat.obj = fixture$seu,
      .meta = sample_meta,
      .sample.var = "sample_id",
      .assay = "RNA",
      .layer.pattern = "counts.",
      .layer.name.source = "missing_column",
      .min.cells.per.sample = 2,
      .verbose = FALSE
    ),
    "not found in Seurat object metadata"
  )
})
