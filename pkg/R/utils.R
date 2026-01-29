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

## ticks ####
### major ticks ####
#' Major tick marks for log scale
#' @export
major_tick <- c(-10000, -1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000, 100000)
names(major_tick) <- major_tick

#' Major tick marks for linear scale
#' @export
major_tick_linear <- c(0, 50000, 100000, 150000, 200000, 250000,20000000)
names(major_tick_linear) <- major_tick_linear # c(0, "50k", "100k", "150k", "200k", "250k")

### minor ticks ####
tick <- NULL
for (m in major_tick) {
  tick <- c(tick, seq(0, 10, by = 1) * m)
}
tick <- unique(tick)

#' Minor tick marks for log scale
#' @export
minor_tick <- sort(tick[!tick %in% major_tick])
minor_tick <- minor_tick[minor_tick >= major_tick[1]]
minor_tick <- minor_tick[minor_tick > -1000]
names(minor_tick) <- minor_tick
rm(tick)

tick <- NULL
for (m in major_tick_linear) {
  tick <- c(tick, seq(0, 10, by = 1) * m)
}
tick <- unique(tick)

#' Minor tick marks for linear scale
#' @export
minor_tick_linear <- sort(tick[!tick %in% major_tick_linear])
minor_tick_linear <- minor_tick_linear[minor_tick_linear >= major_tick_linear[1]]
names(minor_tick_linear) <- minor_tick_linear
rm(tick)

#' Parse channel names or indices from a flow cytometry object
#'
#' @param obj A flow frame object
#' @param channels Channel names or indices
#' @return Numeric indices of the channels
#' @export
# f.parseChannels <- function(obj, channels) {
#   fluos <- parameters(obj)[[1]]
#   markers <- as.vector(parameters(obj)[["desc"]])
#   fluos <- gsub("Comp-", "", fluos)
#   if (missing(channels)) {
#     channels <- c(which(fluos == "SSC-A"), which(fluos == "FSC-A"))
#   } else {
#     channels <- unlist(lapply(channels, function(x) {
#       if (class(x) == "numeric") {
#         idx <- x # using channel index
#       } else if (class(x) == "character") {
#         idx <- c(which(fluos == x), which(markers == x)) #exzact match
#         idx <- idx[duplicated(idx)] # when fluo name equals channel name, choosing either one is the same
#         if (length(idx) == 0) {
#           x <- paste0("\\b", x)
#           idx <- c(grep(x, fluos, ignore.case = T)[1], grep(x, markers, ignore.case = T)[1])
#           idx <- na.omit(idx)
#         } else if (length(idx) > 1) {
#           idx <- idx[endsWith(fluos[idx], "-A")]
#         }
#       }
#       return(idx)
#     }))
#   }
#
#   return(channels)
# }

f.parseChannels <- function (obj, channels)
{
  fluos <- parameters(obj)[[1]]
  markers <- as.vector(parameters(obj)[["desc"]])
  fluos <- gsub("Comp-", "", fluos)

  if (missing(channels)) {
    channels <- c(which(fluos == "SSC-A"), which(fluos ==
                                                   "FSC-A"))
  }
  else {
    channels <- unlist(lapply(channels, function(x) {

      xb <- paste0("\\b", x)
      idx <- c(x, which(fluos == x), which(markers == x), grep(xb, fluos, ignore.case = T)[1],
               grep(xb, markers, ignore.case = T)[1])

      idx <- suppressWarnings(as.numeric(idx))
      idx <- na.omit(idx)

      if(length(idx) >1) {
        if(sum(duplicated(idx))>0) {
          idx <- idx[!duplicated(idx)]
        }
      }
      return(idx)
    }))
  }
  return(channels)
}




#' Format numbers as scientific expressions
#'
#' @param x Numeric values to format
#' @return Formatted scientific expressions
#' @export
scientific <- function(x) { # function for formating the annotation
  ifelse(x %in% c(-1, 1, -10, 10, -100), " ",
         parse(text = gsub(
           "[+]", "",
           gsub(
             "1e", "10^",
             scientific_format()(x)
           )
         )) # skip the -1 -10 1, 10
  )
}

#' Draw axis ticks for flow cytometry plots
#'
#' @param p_chnames Channel names
#' @param major_tick_loc Major tick locations
#' @param minor_tick_loc Minor tick locations
#' @param ranges Plot ranges
#' @param shownumbers Whether to show tick labels
#' @export
f.axis_ticks <- function(p_chnames,
                         major_tick_loc,
                         minor_tick_loc,
                         ranges,
                         shownumbers = FALSE) {
  for (i in 1:2) {
    major.tick <- major_tick_loc[[i]]
    major.tick <- major.tick[major.tick >= ranges[[i]][1] & major.tick <= ranges[[i]][2]]

    if (shownumbers == T) {
      major.tick.labels <- scientific(as.numeric(names(major.tick)))
    } else {
      major.tick.labels <- NA
    }

    axis(i, major.tick,
         col = NA, col.ticks = 1, tck = -0.02,
         labels = major.tick.labels, las = 0, cex.axis = 0.7
    )

    minor.tick <- minor_tick_loc[[i]]
    minor.tick <- minor.tick[minor.tick >= ranges[[i]][1] & minor.tick <= ranges[[i]][2]]
    rug(x = minor.tick, ticksize = -0.01, side = i)
  }
}

#' Add axis labels to flow cytometry plots
#'
#' @param labels Axis labels
#' @param padding Padding for labels
#' @param font Font style
#' @param font.size Font size
#' @export
f.axis_label <- function(labels,
                         padding = 1.3,
                         font = 2,
                         font.size = 0.8) {
  mtext(
    side = 1, line = padding,
    labels[1],
    col = "black",
    adj = 0, # left alighted ,1 center
    font = font,
    cex = font.size
  ) # axis label
  mtext(
    side = 2, line = padding,
    labels[2],
    col = "black", adj = 0,
    font = font,
    cex = font.size
  )
}

#' Transform flow cytometry data
#'
#' @param x Data vector
#' @param method Transformation method
#' @param parameters Transformation parameters
#' @return Transformed data
#' @export
f.transVectData <- function(x,
                            method,
                            parameters = transform_parameters) {
  if (method == "logicle") {
    t <- na.omit(c(parameters["t"], 262143))[1]
    w <- parameters["w"]
    m <- na.omit(c(parameters["m"], 4.5))[1]
    a <- na.omit(c(parameters["a"], 0))[1]
    logicle_trans <- logicleTransform(
      t = t,
      w = parameters["w"],
      m = m,
      a = a,
      "logicle"
    )
    x <- logicle_trans(x)
  } else if (method == "arcsinh") {
    x <- asinh(x / parameters["b"])
  } else if (method == "log") {
    log_trans <- logTransform(logbase = 10, r = 1, d = 1, "log")
    x[x <= 0] <- 0.1

    x <- log_trans(x)
  }

  return(x)
}

#' Transform FSC/SSC data
#'
#' @param x Data vector
#' @param method Transformation method
#' @return Transformed data
#' @export
f.transFscSscData <- function(x, method) {
  if (method == "linear") {
    x
  } else if (method == "log") {
    log_trans <- logTransform(logbase = 10, r = 1, d = 1, "log")
    x[x <= 0] <- 0.1
    x <- log_trans(x)
  }

  return(as.numeric(x))
}

#' Plot gates on flow cytometry plots
#'
#' @param x Gate definitions
#' @param xlim X-axis limits
#' @param ylim Y-axis limits
#' @param log_ssc Whether SSC is log-transformed
#' @export
f.plotGate <- function(x, xlim, ylim, log_ssc, col = "black") {
  lapply(x, function(g) {
    if (class(g)[1] == "matrix") {
      ssc.idx <- which(colnames(g) == "SSC-A")
      if(length(ssc.idx)>0 & log_ssc == T) {
        log_trans <- logTransform(logbase = 10, r = 1, d = 1, "log")
        g[g[,ssc.idx]<= 0,] <- 0.1
        g[,ssc.idx] <- log_trans( g[,ssc.idx])
      }

      lines(g,
            type = "l",
            lwd = 2,
            col = col
      )
    } else if (class(g)[1] == "numeric") {
      abline(
        v = g[1], h = g[2],
        lty = 1,
        lwd = 2,
        col = col
      )
    }
  })
}

#' Map numeric values to color gradients
#'
#' @param numbers Numeric values
#' @param palette Color palette
#' @param lower_percentile Lower percentile for normalization
#' @param upper_percentile Upper percentile for normalization
#' @param breaks Number of color bins
#' @return Character vector of hex color codes
#' @export
number_to_color <- function(numbers, palette = "magma", lower_percentile = 0.01, upper_percentile = 0.99, breaks = 256) {
  # Validate the palette input
  if (!palette %in% c("magma", "inferno", "plasma", "viridis", "cividis")) {
    stop("Invalid palette. Choose from 'magma', 'inferno', 'plasma', 'viridis', or 'cividis'.")
  }

  # Compute robust lower and upper bounds using percentiles
  lower_bound <- quantile(numbers, probs = lower_percentile, na.rm = TRUE)
  upper_bound <- quantile(numbers, probs = upper_percentile, na.rm = TRUE)

  # Normalize the numbers to the range [0, 1] using robust bounds
  normalized_numbers <- (numbers - lower_bound) / (upper_bound - lower_bound)

  # Clip values outside the range [0, 1] to avoid extrapolation
  normalized_numbers <- pmin(pmax(normalized_numbers, 0), 1)

  # Map the normalized numbers to the selected color palette
  color_function <- switch(palette,
                           magma = viridis::magma,
                           inferno = viridis::inferno,
                           plasma = viridis::plasma,
                           viridis = viridis::viridis,
                           cividis = viridis::cividis)

  colors <- color_function(breaks)[as.numeric(cut(normalized_numbers, breaks = breaks))]

  return(colors)
}




#' Color palette with 103 distinct colors
#' @export
Color103 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
  "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
  "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
  "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
  "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
  "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
  "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
  "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
  "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
  "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
  "#1D1702", "#365D25"
)
names(Color103) <- seq(length(Color103))

