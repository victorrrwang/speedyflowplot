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

#' Create a 2D density plot for flow cytometry data
#'
#' A simple implementation of scattermore functionality for flow cytometry data visualization.
#'
#' @param obj A flow frame, matrix, or data frame containing flow cytometry data
#' @param channels Channels to plot (column indices or names)
#' @param title Plot title
#' @param title.font.size Font size for the title
#' @param res Resolution of the plot
#' @param color Color for the points
#' @param xlim X-axis limits
#' @param ylim Y-axis limits
#' @param plot_ceiling Maximum value for plotting
#' @param density.overlay Whether to overlay density plots
#' @param dens.color Color for density overlay
#' @param dens.alpha Alpha transparency for density overlay
#'
#' @return A plot of the flow cytometry data
#' @export
flowPlot <- function(obj,
                     channels,
                     title = "",
                     title.font.size = 1.2,
                     res = 256,
                     color = "#00000010",
                     xlim = NA,
                     ylim = NA,
                     plot_ceiling = NA,
                     density.overlay = TRUE, 
                     dens.color = "red",
                     dens.alpha = 0.5,
                     show.tick.labels = TRUE
                     ) {
  if (class(obj)[1] == "matrix") {
    dat <- obj[, channels]
  } else if (class(obj)[1] == "data.frame") {
    dat <- as.matrix(obj[, channels])
  } else if (class(obj)[1] == "flowFrame") {
    dat <- exprs(obj)[, channels]
  } else if (class(obj)[1] == "CellPopulation") {
    dat <- exprs(obj@flow.frame[obj@index, channels])
  }
  
  if(nrow(dat) > 0) {
    if(nrow(dat) > 500000) {
      set.seed(23)
      dat <- dat[sample(nrow(dat), 500000), ]
    }
    
    # set limits
    if (is.na(plot_ceiling)) {
      plot_ceiling <- max(dat)
    }
    
    if (is.na(sum(xlim))) {
      xlim <- range(dat[, 1])
    } else {
      dat[, 1][dat[, 1] < xlim[1]] <- xlim[1]
      dat[, 1][dat[, 1] > xlim[2]] <- xlim[2]
    }
    if (is.na(sum(ylim))) {
      ylim <- range(dat[, 2])
    } else {
      dat[, 2][dat[, 2] < ylim[1]] <- ylim[1]
      dat[, 2][dat[, 2] > ylim[2]] <- ylim[2]
    }
    
    if(missing(res)) {
      cn <- nrow(dat)
      
      if(cn >= 5000) {
        res = 256
      } else if (cn < 5000 & cn >= 2500) {
        res = 196
      } else if(cn < 2500 & cn >= 200){
        res = 128
      } else {
        res = 64
      }
    }
    
    par(pty = "s") # make the plot square
    
    rgbwt <- scatter_points_rgbwt(dat, 
                                  xlim = xlim,
                                  ylim = ylim,
                                  out_size = c(res, res), 
                                  RGBA = col2rgb(color, alpha=TRUE))
    rgbwt[,,5] <- 1-rgbwt[,,5]
    rgbwt[,,5][rgbwt[,,5] == 0] <- 1
    rstr <- rgba_int_to_raster(rgbwt_to_rgba_int(rgbwt))
    
    plot(c(),
         xlim = xlim, xlab = channels[1], 
         ylim = ylim, ylab = channels[2], 
         xaxt = if(show.tick.labels) "s" else "n",  # Add this
         yaxt = if(show.tick.labels) "s" else "n"    # Add this
    )
    title(title, adj = 0, line = 0.3, cex.main = title.font.size)
    rasterImage(rstr,
                xleft = xlim[1],
                xright = xlim[2],
                ybottom = ylim[1],
                ytop = ylim[2], 
                interpolate = FALSE
    )
    
    if (density.overlay & nrow(dat) > 10) {
      adjust.dens = 1
      x.dens <- density(dat[, 1], adjust = adjust.dens)
      y.dens <- density(dat[, 2], adjust = adjust.dens)
      
      x.axis <- x.dens$x
      y.axis <- y.dens$x
      x.dens <- scales::rescale(x.dens$y, ylim)
      y.dens <- scales::rescale(y.dens$y, xlim)
      
      lines(x.axis, x.dens, col = scales::alpha(dens.color, dens.alpha))
      lines(y.dens, y.axis, col = scales::alpha(dens.color, dens.alpha))
    }
  } else {
    plot(c(),
         xlim = c(0,4.5), xlab = channels[1], 
         ylim = c(0,4.5), ylab = channels[2], 
    )
  }
  
}
