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

#' Fetch trajectory simulation dataset from miloR
#'
#' Downloads and extracts the sim_trajectory dataset from the miloR package GitHub repository.
#' This dataset contains simulated single-cell trajectory data with two conditions (A and B)
#' across three replicates each.
#'
#' @return A list containing:
#'   \item{meta}{Sample metadata with Condition, Replicate, and Sample columns}
#'   \item{SCE}{SingleCellExperiment object containing the count data}
#'
#' @examples
#' \dontrun{
#' # Fetch the trajectory data
#' sim_trajectory_data <- fetch_trajectory_data()
#' 
#' # Extract components
#' sim_trajectory.meta <- sim_trajectory_data$meta
#' sim_trajectory <- sim_trajectory_data$SCE
#' }
#'
#' @source 
#' Data originally from: Dann, E., Henderson, N.C., Teichmann, S.A. et al. 
#' Differential abundance testing on single-cell data using k-nearest neighbor graphs. 
#' Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01033-z
#' 
#' miloR GitHub repository: https://github.com/MarioniLab/miloR
#' Direct data link: https://github.com/MarioniLab/miloR/blob/bdecaebeb5595545f9b9c8f2defb321519d98a70/data/sim_trajectory.RData
#'
#' @note The downloaded data is subject to GPL v3 license terms from miloR package
#'
#' @export
fetch_trajectory_data <- function() {
  
  # Check if we have internet connection
  if (!curl::has_internet()) {
    stop("Internet connection required to fetch trajectory data from miloR repository")
  }
  
  # URL to the data file
  data_url <- "https://github.com/MarioniLab/miloR/raw/bdecaebeb5595545f9b9c8f2defb321519d98a70/data/sim_trajectory.RData"
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".RData")
  
  # Download the data
  message("Downloading trajectory data from miloR repository...")
  utils::download.file(url = data_url, destfile = temp_file, mode = "wb", quiet = TRUE)
  
  # Load the data
  env <- new.env()
  load(temp_file, envir = env)
  
  # Clean up
  unlink(temp_file)
  
  # Extract the sim_trajectory object
  sim_trajectory <- env$sim_trajectory
  
  # add meta data to SCE
  SummarizedExperiment::colData(x = sim_trajectory$SCE) <-
    as.list(x = sim_trajectory$meta) |>
    S4Vectors::DataFrame()
  
  colnames(x = sim_trajectory$SCE) <-
    sim_trajectory$meta$cell_id
  
  message("Successfully fetched trajectory data with ", ncol(sim_trajectory$SCE), " cells and ", nrow(sim_trajectory$SCE), " features")
  
  return(sim_trajectory)
}
