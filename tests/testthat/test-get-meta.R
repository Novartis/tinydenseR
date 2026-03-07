# ======================================================================
# get.meta.HDF5AnnData tests
# ======================================================================

test_that("get.meta.HDF5AnnData extracts sample-level metadata", {
  skip_on_cran()
  skip_if_not_installed("anndataR")
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)

  adata <- anndataR::as_AnnData(x = sce, output_class = "HDF5AnnData",
                                file = h5ad_path)

  meta <- get.meta.HDF5AnnData(.h5ad.obj = adata,
                                .sample.var = "Sample",
                                .verbose = FALSE)

  expect_s3_class(meta, "data.frame")

  n_samples <- length(unique(as.character(adata$obs[["Sample"]])))
  expect_equal(nrow(meta), n_samples)

  expect_true(all(rownames(meta) %in% as.character(adata$obs[["Sample"]])))
  expect_true("Sample" %in% colnames(meta))
})

test_that("get.meta.HDF5AnnData excludes cell-level columns", {
  skip_on_cran()
  skip_if_not_installed("anndataR")
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)

  adata <- anndataR::as_AnnData(x = sce, output_class = "HDF5AnnData",
                                file = h5ad_path)

  meta <- get.meta.HDF5AnnData(.h5ad.obj = adata,
                                .sample.var = "Sample",
                                .verbose = FALSE)

  obs_df <- as.data.frame(adata$obs)
  for (col in colnames(meta)) {
    vals_per_sample <- tapply(obs_df[[col]], obs_df[["Sample"]],
                              function(x) length(unique(x)))
    expect_true(all(vals_per_sample == 1),
                info = paste("Column", col, "is not sample-level"))
  }
})

test_that("get.meta.HDF5AnnData errors on invalid inputs", {
  skip_if_not_installed("anndataR")

  expect_error(get.meta.HDF5AnnData(.h5ad.obj = data.frame(),
                                     .sample.var = "x"),
               "must be.*HDF5AnnData")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)

  adata <- anndataR::as_AnnData(x = sce, output_class = "HDF5AnnData",
                                file = h5ad_path)

  expect_error(get.meta.HDF5AnnData(.h5ad.obj = adata,
                                     .sample.var = 123),
               "must be a single character string")
})

test_that("get.meta dispatches to HDF5AnnData", {
  skip_on_cran()
  skip_if_not_installed("anndataR")
  skip_if_not(curl::has_internet(), "No internet connection")

  trajectory_data <- fetch_trajectory_data()
  sce <- trajectory_data$SCE

  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)

  adata <- anndataR::as_AnnData(x = sce, output_class = "HDF5AnnData",
                                file = h5ad_path)

  meta <- get.meta(.obj = adata, .sample.var = "Sample",
                   .verbose = FALSE)

  expect_s3_class(meta, "data.frame")
  n_samples <- length(unique(as.character(adata$obs[["Sample"]])))
  expect_equal(nrow(meta), n_samples)
})
