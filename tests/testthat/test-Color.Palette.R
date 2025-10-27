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

test_that("Color.Palette is a matrix with correct dimensions and names", {
  expect_true(exists("Color.Palette"))
  expect_true(is.matrix(Color.Palette))
  expect_equal(dim(Color.Palette), c(6, 7))
  expect_true(!is.null(dimnames(Color.Palette)))
})
