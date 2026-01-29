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

#' Blend Multiple RGBA Arrays with Selective Transparency
#'
#' @description
#' Blends multiple RGBA (Red, Green, Blue, Alpha) arrays with a special transparency
#' rule where only overlapping pixels have reduced transparency. This creates a
#' visual effect where overlapping points are more transparent than non-overlapping ones.
#'
#' @param fRGBA_list A list of RGBA arrays. Each array should have dimensions
#'   [height, width, 4] where the fourth dimension represents RGBA channels.
#'   Values should be in float format (0-1 range).
#'
#' @return An array with dimensions [height, width, 4] containing the blended RGBA values.
#'
#' @details
#' The function performs alpha blending with these key features:
#' \itemize{
#'   \item Identifies overlapping pixels (where both source and destination alpha > 0)
#'   \item Reduces transparency by 80% only for overlapping pixels
#'   \item Preserves original transparency for non-overlapping pixels
#'   \item Blends colors according to modified alpha values
#' }
#'
#' The blending formula used is:
#' \deqn{C_{out} = C_1 \alpha_1 + C_2(1-\alpha_1)}
#' where \eqn{C} represents color channels (R,G,B) and \eqn{\alpha} represents opacity.
#'
#' @examples
#' # Create two sample RGBA arrays
#' dims <- c(100, 100, 4)
#' rgba1 <- array(0, dims)
#' rgba2 <- array(0, dims)
#'
#' # Set some test values
#' rgba1[40:60, 40:60, ] <- c(1, 0, 0, 0.5)  # Red square
#' rgba2[50:70, 50:70, ] <- c(0, 0, 1, 0.5)  # Blue square
#'
#' # Blend the arrays
#' result <- blend_rgba_float3(list(rgba1, rgba2))
#'
#' @export
#' @useDynLib speedyflowplot, .registration=TRUE
blend_rgba_float3 <- function (fRGBA_list, alpha = 0.8) 
{
    if (length(fRGBA_list) < 1) 
        stop("No input RGBA given.")
    if (length(fRGBA_list) == 1) 
        return(fRGBA_list[[1]])
    fRGBA_1 <- fRGBA_list[[1]]
    if (!is.array(fRGBA_1) || dim(fRGBA_1)[3] != 4) 
        stop("unsupported RGBA format")
    for (i in 2:length(fRGBA_list)) {
        fRGBA_2 <- fRGBA_list[[i]]
        if (!is.array(fRGBA_2) || dim(fRGBA_2)[3] != 4) 
            stop("unsupported RGBA format")
        if ((dim(fRGBA_1)[1] != dim(fRGBA_2)[1]) || (dim(fRGBA_1)[2] != 
            dim(fRGBA_2)[2])) 
            stop("input bitmap dimensions differ")
        rows <- dim(fRGBA_1)[1]
        cols <- dim(fRGBA_1)[2]
        overlap_mask <- (fRGBA_1[, , 4] > 0 & fRGBA_2[, , 4] > 
            0)
        A_1 <- fRGBA_1[, , 4]
        A_1[overlap_mask] <- A_1[overlap_mask] * alpha
        A_2 <- fRGBA_2[, , 4]
        fRGBA <- array(0, c(rows, cols, 4))
        fRGBA[, , 1] <- fRGBA_1[, , 1] * A_1 + fRGBA_2[, , 1] * 
            (1 - A_1)
        fRGBA[, , 2] <- fRGBA_1[, , 2] * A_1 + fRGBA_2[, , 2] * 
            (1 - A_1)
        fRGBA[, , 3] <- fRGBA_1[, , 3] * A_1 + fRGBA_2[, , 3] * 
            (1 - A_1)
        fRGBA[, , 4] <- A_1 + (A_2 * (1 - A_1))
        fRGBA_1 <- fRGBA
    }
    return(fRGBA_1)
}