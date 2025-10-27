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

# Test for leiden.cluster

test_that("leiden.cluster returns integer cluster assignments", {
  .sim.matrix <- matrix(data = runif(n = 400), nrow = 20, ncol = 20)
  .sim.matrix <- (.sim.matrix + t(.sim.matrix)) / 2 # make symmetric
  result <- tinydenseR:::leiden.cluster(.sim.matrix = .sim.matrix,
                          .resolution.parameter = 1.0)
  expect_type(object = result, type = "integer")
  expect_length(object = result, n = 20)
})

# Test for get.adj.matrix

test_that("get.adj.matrix returns a sparse matrix with correct dimensions", {
  .nn.idx <- matrix(data = sample(x = 1:10, size = 30, replace = TRUE), nrow = 10, ncol = 3)
  result <-  tinydenseR:::get.adj.matrix(.nn.idx = .nn.idx)
  expect_s4_class(object = result, class = "dgCMatrix")
  expect_equal(object = dim(x = result), expected = c(10, 10))
})

# Test for fast.jaccard.r

test_that("fast.jaccard.r returns a symmetric sparse matrix", {
  .adj.matrix <- Matrix::rsparsematrix(nrow = 10, ncol = 10, density = 0.2)
  result <-  tinydenseR:::fast.jaccard.r(.adj.matrix = .adj.matrix,
                           .prune = 1/15)
  expect_s4_class(object = result, class = "generalMatrix")
  expect_equal(object = dim(x = result), expected = c(10, 10))
  expect_true(Matrix::isSymmetric(result))
})

# Test for get.graph

test_that("get.graph adds graph to .lm.obj", {
  .cells <- list(sample1 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample1",
                                                         formatC(x = 1:90,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8"))),
                 sample2 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample2",
                                                         formatC(x = 1:90,
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
  result <- get.landmarks(.lm.obj = result,
                           .verbose = FALSE)
  result <- get.graph(.lm.obj = result,
                        .k = 3,
                        .scale = FALSE,
                        .verbose = FALSE,
                        .seed = 42)
  expect_type(object = result, type = "list")
  expect_true("graph" %in% names(x = result))
  expect_true("uwot" %in% names(x = result$graph))
})

# Test for get.map

test_that("get.map adds mapping info to .lm.obj", {
  .cells <- list(sample1 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample1",
                                                         formatC(x = 1:90,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8"))),
                 sample2 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample2",
                                                         formatC(x = 1:90,
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
  result <- get.landmarks(.lm.obj = result,
                          .verbose = FALSE)
  result <- get.graph(.lm.obj = result,
                         .k = 3,
                         .scale = FALSE,
                         .verbose = FALSE,
                         .seed = 42)
  expect_type(object = result, type = "list")
})

# Test for get.lm.features.stats

test_that("get.lm.features.stats returns a list of feature stats", {
  .cells <- list(sample1 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample1",
                                                         formatC(x = 1:90,
                                                                 width = 2,
                                                                 format = "d",
                                                                 flag = "0")), 
                                                  c("CD3", "CD4", "CD8"))),
                 sample2 = matrix(data = runif(n = 270), 
                                  nrow = 90, 
                                  ncol = 3, 
                                  dimnames = list(paste0("sample2",
                                                         formatC(x = 1:90,
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
  result <- get.landmarks(.lm.obj = result,
                          .verbose = FALSE)
  result <- get.graph(.lm.obj = result,
                         .k = 3,
                         .scale = FALSE,
                         .verbose = FALSE,
                         .seed = 42)
  expect_type(object = result, type = "list")
  result <- get.lm.features.stats(.lm.obj = result)
  expect_type(object = result, type = "list")
})
