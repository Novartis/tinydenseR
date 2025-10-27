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

#' Color Palette
#'
#' This is a color palette
#'
#' @format A matrix with 6 rows and 7 columns
#' @export
#' @name Color.Palette

# Color for emphasis or distinction

Color.Palette <- structure(
  c(  "#0460A9", "#CDDFEE", "#9BBFDD", "#68A0CB", "#03487F", "#023054",
      "#E74A21", "#FADBD3", "#F5B7A6", "#F1927A", "#AD3819", "#742510",
      "#EC9A1E", "#FBEBD2", "#F7D7A5", "#F4C278", "#B17416", "#764D0F",
      "#8D1F1B", "#E8D2D1", "#D1A5A4", "#BB7976", "#6A1714", "#46100E",
      "#7F7F7F", "#E5E5E5", "#CCCCCC", "#B2B2B2", "#5F5F5F", "#404040",
      "#CCCCCC", "#F5F5F5", "#EBEBEB", "#E0E0E0", "#999999", "#666666",
      "#404040", "#D9D9D9", "#B3B3B3", "#8C8C8C", "#303030", "#202020"
  ),
  .Dim = c(6L, 7L),
  .Dimnames = list(c("Hue", "Tint3", "Tint2", "Tint1", "Shade1", "Shade2"),
                   c("Blue", "Sienna", "Apricot", "Carmine", "Gray", "Light Gray", "Dark Gray")))
