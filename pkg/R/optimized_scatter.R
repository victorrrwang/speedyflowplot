# This file is part of speedyflowplot.
#
# Copyright (C) 2024 Victor Wang <VictorrWang@gmail.com>
#
# speedyflowplot is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# speedyflowplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with speedyflowplot. If not, see <https://www.gnu.org/licenses/>.

#' scatter_points_optimized
#'
#' Optimized version of scatter_points_rgbwt that uses raw uint8_t arrays for faster performance.
#'
#' @param xy 2-column matrix with N point coordinates (X and Y) in rows.
#'
#' @param xlim,ylim 2-element vector of rendered area limits (position of the first pixel on the
#'                   left/top, and the last pixel on the right/bottom).
#'                   You can flip the image coordinate system by flipping the `*lim` vectors.
#'
#' @param out_size 2-element vector size of the result raster, defaults to `c(512L,512L)`.
#'
#' @param RGBA Point colors. Either a 4-element vector that specifies the same color for all points,
#'             or 4-by-N matrix that specifies color for each of the individual points.
#'             Color is specified using integer RGBA; i.e. the default black is `c(0,0,0,255)`.
#'
#' @param map Vector with N integer indices to `palette`. Overrides RGBA-based coloring.
#'
#' @param palette Matrix 4-by-K matrix of RGBA colors used as a palette lookup for the `map`
#'                that gives the point colors. K is at least `max(map)`.
#'                Notably, using a palette may be faster than filling and processing the whole RGBA matrix.
#'
#' @return A raw vector containing RGBA bitmap data.
#'
#' @export
#' @useDynLib speedyflowplot, .registration=TRUE
scatter_points_optimized <- function(xy,
                                    xlim = c(min(xy[, 1]), max(xy[, 1])),
                                    ylim = c(min(xy[, 2]), max(xy[, 2])),
                                    out_size = c(512, 512),
                                    RGBA = c(0, 0, 0, 255),
                                    map = NULL,
                                    palette = NULL) {
  if (!is.numeric(xlim) || length(xlim) != 2) stop("invalid xlim")
  if (!is.numeric(ylim) || length(ylim) != 2) stop("invalid ylim")
  if (!is.numeric(out_size) || length(out_size) != 2) stop("invalid out_size")

  if (dim(xy)[2] != 2) stop("2-column xy input expected")
  n <- dim(xy)[1]

  size_x <- as.integer(out_size[1])
  size_y <- as.integer(out_size[2])
  
  # Create raw output buffer (RGBA for each pixel)
  raw_output <- raw(size_x * size_y * 4)
  
    # Call the appropriate C function based on the input
    if (is.numeric(map)) {
      if (length(map) != n) stop("wrong size of map")
      if (any(map < 1)) stop("indices in map must start from 1")
      if (!is.matrix(palette) && !is.array(palette)) stop("unsupported palette format")
      if (dim(palette)[1] != 4) stop("unsupported palette format")
      if (max(map) > dim(palette)[2]) stop("map indices too high for this palette")
      
      # Use the existing .C interface instead of .Call
      # Create a temporary RGBWT array
      rgbwt <- array(0, c(size_y, size_x, 5))
      rgbwt[, , 5] <- 1 # initialize transparency (multiplying)
      
      # Call the indexed version using .C
      result <- .C("scatter_indexed_rgbwt",
                  dimen = as.integer(c(size_x, size_y, n)),
                  xlim = as.single(xlim),
                  ylim = as.single(ylim),
                  palette = as.single(palette / 255),
                  fRGBWT = as.single(rgbwt),
                  map = as.integer(map - 1L), # Convert to 0-based indexing
                  xy = as.single(xy))
      
      # Convert RGBWT to RGBA
      rgbwt <- array(result$fRGBWT, c(size_y, size_x, 5))
      rgba <- rgbwt_to_rgba_int(rgbwt)
      
      # Copy to raw output
      for (i in 1:4) {
        raw_output[seq(i, length(raw_output), by = 4)] <- as.raw(rgba[,,i])
      }
    } else {
      # Prepare RGBA data
      if (is.vector(RGBA) || ((is.matrix(RGBA) || is.array(RGBA)) && dim(RGBA)[2] == 1)) {
        if (length(RGBA) != 4) stop("RGBA vector of length 4 expected")
        rgba_data <- RGBA
      } else if (is.matrix(RGBA) || is.array(RGBA)) {
        if (dim(RGBA)[1] != 4) stop("RGBA matrix with 4 rows expected")
        if (dim(RGBA)[2] != n) stop("incorrect number of colors in RGBA")
        rgba_data <- RGBA
      } else {
        stop("unsupported coloring type")
      }
      
      # Use the existing .C interface instead of .Call
      # Create a temporary RGBWT array
      rgbwt <- array(0, c(size_y, size_x, 5))
      rgbwt[, , 5] <- 1 # initialize transparency (multiplying)
      
      # Call the single color version using .C
      if (is.vector(rgba_data) || ((is.matrix(rgba_data) || is.array(rgba_data)) && dim(rgba_data)[2] == 1)) {
        result <- .C("scatter_singlecolor_rgbwt",
                    dimen = as.integer(c(size_x, size_y, n)),
                    xlim = as.single(xlim),
                    ylim = as.single(ylim),
                    RGBA = as.single(rgba_data / 255),
                    fRGBWT = as.single(rgbwt),
                    xy = as.single(xy))
      } else {
        # Call the multi-color version using .C
        result <- .C("scatter_multicolor_rgbwt",
                    dimen = as.integer(c(size_x, size_y, n)),
                    xlim = as.single(xlim),
                    ylim = as.single(ylim),
                    RGBA = as.single(rgba_data / 255),
                    fRGBWT = as.single(rgbwt),
                    xy = as.single(xy))
      }
      
      # Convert RGBWT to RGBA
      rgbwt <- array(result$fRGBWT, c(size_y, size_x, 5))
      rgba <- rgbwt_to_rgba_int(rgbwt)
      
      # Copy to raw output
      for (i in 1:4) {
        raw_output[seq(i, length(raw_output), by = 4)] <- as.raw(rgba[,,i])
      }
    }
  
  return(raw_output)
}

#' apply_kernel_optimized
#'
#' Optimized version of apply_kernel_rgbwt that uses raw uint8_t arrays for faster performance.
#'
#' @param raw_input Raw vector containing RGBA bitmap data.
#'
#' @param out_size 2-element vector size of the bitmap.
#'
#' @param filter Use the pre-defined filter, either `circle`, `square`, `gauss`. Defaults to `circle`.
#'
#' @param mask Custom kernel used for blurring, overrides `filter`. Must be a square matrix of odd size.
#'
#' @param radius Radius of the kernel (counted without the "middle" pixel"), defaults to 2. The generated kernel matrix will be a square with (2*radius+1) pixels on each side.
#'
#' @param sigma Radius of the Gaussian function selected by `filter`, defaults to `radius/2`.
#'
#' @return A raw vector containing the blurred RGBA bitmap data.
#'
#' @export
#' @useDynLib speedyflowplot, .registration=TRUE
apply_kernel_optimized <- function(raw_input,
                                  out_size,
                                  filter = "circle",
                                  mask = default_kernel(filter, radius, sigma),
                                  radius = 2,
                                  sigma = radius / 2) {
  if (!is.raw(raw_input)) stop("raw_input must be a raw vector")
  if (!is.numeric(out_size) || length(out_size) != 2) stop("invalid out_size")
  
  size_x <- as.integer(out_size[1])
  size_y <- as.integer(out_size[2])
  
  if (length(raw_input) != size_x * size_y * 4) stop("raw_input size doesn't match out_size")
  
  if (!is.matrix(mask) && !is.array(mask)) stop("kernel in matrix or array form expected")
  if (dim(mask)[1] != dim(mask)[2]) stop("kernel in square matrix expected")
  if (dim(mask)[1] %% 2 == 0) stop("kernel with odd size expected")
  
  # Convert raw input to RGBWT
  # First convert raw to RGBA
  rgba_matrix <- matrix(as.integer(raw_input), nrow = 4)
  rgba_array <- array(0, c(size_y, size_x, 4))
  for (i in 1:4) {
    rgba_array[,,i] <- matrix(rgba_matrix[i,], nrow = size_y, ncol = size_x, byrow = TRUE)
  }
  
  # Convert RGBA to RGBWT
  rgbwt <- array(0, c(size_y, size_x, 5))
  rgbwt[,,1:3] <- rgba_array[,,1:3] * rgba_array[,,4] / 255
  rgbwt[,,4] <- rgba_array[,,4] / 255
  rgbwt[,,5] <- 1 - rgbwt[,,4]
  
  # Create blurred RGBWT array
  blurred_rgbwt <- array(0, c(size_y, size_x, 5))
  blurred_rgbwt[,,5] <- 1 # initialize transparency (multiplying)
  
  # Call the kernel function using .C
  kernel_pixels <- floor(dim(mask)[1] / 2)
  result <- .C("kernel_rgbwt",
              dimen = as.integer(c(size_x, size_y, kernel_pixels, 0)), # 0 threads means auto-detect
              kernel = as.single(mask),
              blurred_fRGBWT = as.single(blurred_rgbwt),
              fRGBWT = as.single(rgbwt))
  
  # Convert blurred RGBWT to RGBA
  blurred_rgbwt <- array(result$blurred_fRGBWT, c(size_y, size_x, 5))
  rgba <- rgbwt_to_rgba_int(blurred_rgbwt)
  
  # Copy to raw output
  raw_output <- raw(size_x * size_y * 4)
  for (i in 1:4) {
    raw_output[seq(i, length(raw_output), by = 4)] <- as.raw(rgba[,,i])
  }
  
  return(raw_output)
}

#' raw_to_raster
#'
#' Convert a raw RGBA bitmap to a raster object for plotting.
#'
#' @param raw_data Raw vector containing RGBA bitmap data.
#'
#' @param out_size 2-element vector size of the bitmap.
#'
#' @return A raster object.
#'
#' @export
raw_to_raster <- function(raw_data, out_size) {
  if (!is.raw(raw_data)) stop("raw_data must be a raw vector")
  if (!is.numeric(out_size) || length(out_size) != 2) stop("invalid out_size")
  
  size_x <- as.integer(out_size[1])
  size_y <- as.integer(out_size[2])
  
  if (length(raw_data) != size_x * size_y * 4) stop("raw_data size doesn't match out_size")
  
  # Convert raw data to a matrix of colors
  rgba_matrix <- matrix(as.integer(raw_data), nrow = 4)
  
  # Create color strings
  colors <- rgb(rgba_matrix[1,], rgba_matrix[2,], rgba_matrix[3,], rgba_matrix[4,], maxColorValue = 255)
  
  # Create raster - fix the orientation issue
  color_matrix <- matrix(colors, nrow = size_y, ncol = size_x)
  # Flip the matrix to match the original orientation
  raster <- color_matrix[size_y:1, ]
  
  return(as.raster(raster))
}
