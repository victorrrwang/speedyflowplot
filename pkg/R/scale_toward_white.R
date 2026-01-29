# This file is part of speedyflowplot.
#
# Copyright (C) 2022 Mirek Kratochvil <exa.exa@gmail.com>
#               2022-2023 Tereza Kulichova <kulichova.t@gmail.com>
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

#' Scale RGB Colors Toward White Based on Density
#'
#' @param rgbwt An array of dimensions [height, width, 5] containing RGBWT values. 
#'   The first three channels are RGB colors, the fourth is weight/density, and the fifth is transparency.
#' @param base_color A numeric vector of length 3 specifying the RGB base color (values between 0-255) 
#'   from which the scaling toward white begins.
#' @param whiteness A numeric value between 0 and 1 specifying the intensity of the whitening effect. 
#'   Default is 0.9. Higher values create more intense whitening.
#' @param gamma Numeric value for gamma correction. If NA (default), logarithmic scaling is used. 
#'   If a number is provided, it applies non-linear scaling with the specified gamma value.
#'
#' @return Returns an array of the same dimensions as the input rgbwt, with modified RGB values 
#'   scaled toward white based on the density channel. The alpha channel is set to 255 (fully opaque).
#'
#' @examples
#' # Create sample RGBWT array
#' height <- 100
#' width <- 100
#' rgbwt <- array(0, dim = c(height, width, 5))
#' rgbwt[40:60, 40:60, 4] <- 1  # Set some weights in the center
#'
#' # Apply whitening effect with blue base color
#' base_color <- c(0, 0, 255)  # Blue
#' result <- scale_toward_white(rgbwt, base_color)
#'
#' @export
#' @useDynLib speedyflowplot, .registration=TRUE
scale_toward_white <- function (rgbwt, base_color, whiteness = 0.9, gamma = NA) 
{
    whiteness <- pmin(1, pmax(0, whiteness))
    white <- c(255, 255, 255)
    whitened <- round(base_color + whiteness * (white - base_color))
    white_color <- pmin(255, pmax(0, whitened))
    weight <- rgbwt[, , 4]
    max_weight <- max(weight, na.rm = TRUE)
    if (is.na(gamma)) {
        density <- log(weight + 1)/log(max_weight + 1)
    }
    else {
        density <- (weight/max_weight)^gamma
    }
    new_rgbwt <- rgbwt
    for (i in 1:3) {
        new_rgbwt[, , i] <- pmin(255, round(base_color[i] + density * 
            (white_color[i] - base_color[i])))
    }
    new_rgbwt[, , 4] <- 255
    return(new_rgbwt)
}