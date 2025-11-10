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

# Test for setup.lm.obj

test_that("setup.lm.obj returns a list with correct names", {
  .cells <- list(sample1 = matrix(data = runif(n = 30), 
                                  nrow = 10, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample1",
                                                         formatC(x = 1:10,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8"))),
                 sample2 = matrix(data = runif(n = 30), 
                                  nrow = 10, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample2",
                                                         formatC(x = 1:10,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8")))) |>
    lapply(FUN = function(x){
      uri <- tempfile(fileext = ".RDS")
      
      saveRDS(object = x,
              file = uri,
              compress = FALSE)
      
      return(uri)
    })
  .meta <- data.frame(row.names = c("sample1", "sample2"),
                      group = c("A", "B"))
  result <- setup.lm.obj(.cells = .cells,
                         .meta = .meta,
                         .markers = c("CD3", "CD4", "CD8"),
                         .assay.type = "cyto",
                         .verbose = FALSE)
  expect_type(object = result, type = "list")
  expect_true(all(c("cells",
                    "lm",
                    "scaled.lm",
                    "raw.lm",
                    "metadata",
                    "key",
                    "pca",
                    "graph",
                    "map",
                    "assay.type",
                    "spec",
                    "n.threads",
                    "markers",
                    "harmony.var",
                    "interact.plot",
                    "harmony.obj") %in% names(x = result)))
})

# Test error conditions for setup.lm.obj

test_that("setup.lm.obj throws error when .cells names don't match .meta rownames", {
  .cells <- list(wrong_name = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.lm.obj(.cells = .cells, .meta = .meta),
               "Sample names mismatch between .cells and .meta")
})

test_that("setup.lm.obj throws error when .harmony.var not in metadata", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.lm.obj(.cells = .cells, .meta = .meta, .harmony.var = "missing_var"),
               "Variables not found in metadata")
})

test_that("setup.lm.obj throws error when .markers insufficient for cyto", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.lm.obj(.cells = .cells, .meta = .meta, .markers = c("CD3"), .assay.type = "cyto"),
               ".markers must contain at least 3 markers for meaningful dimensionality reduction")
})

test_that("setup.lm.obj throws error when .markers used with RNA", {
  .cells <- list(sample1 = tempfile())
  .meta <- data.frame(row.names = "sample1", group = "A")
  expect_error(setup.lm.obj(.cells = .cells, .meta = .meta, .markers = c("CD3"), .assay.type = "RNA"),
               ".markers argument only applies to cytometry data")
})

# Test for get.landmarks

test_that("get.landmarks validates input properly", {
  # Test that get.landmarks expects a proper .lm.obj structure
  expect_error(get.landmarks(.lm.obj = list()),
               "Invalid .lm.obj: missing structure")  # Should error due to missing required fields
  
  # Test with incomplete object
  incomplete_obj <- list(
    assay.type = "RNA",
    cells = list()
  )
  expect_error(get.landmarks(.lm.obj = incomplete_obj),
               "Invalid .lm.obj structure")  # Should error due to missing required fields
})
