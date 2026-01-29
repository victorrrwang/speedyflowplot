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

#' Advanced flow cytometry plotting function
#'
#' Create advanced visualizations of flow cytometry data with multiple options for transformations,
#' population coloring, and density overlays.
#'
#'
#'
# Import the transformList class from flowCore
#' @importClassesFrom flowCore transformList
#'
#' # Define the method with proper signature
setMethod("summary",
          signature=signature(object="transformList"),
          definition=function(object,...)
          {
            sapply(object@transforms,function(transMap)
              summary(new("transform", .Data =
                            transMap@f))
              ,simplify =FALSE)
          })


#'
#'
#' @param obj A flow frame, GatingSet, or FlowSet object
#' @param channels Channels to be plotted
#' @param colormap_marker Channel to use for color mapping
#' @param pop_idex.list List of cell indices for different populations
#' @param plotting_background Whether to plot background data
#' @param title Plot title
#' @param title.font.size Font size for title
#' @param comp Whether to apply compensation
#' @param spillovermatix Spillover matrix for compensation
#' @param min.xy Minimum values for x and y axes
#' @param max.xy Maximum values for x and y axes
#' @param lowerlimit.quantile Lower limit quantile for trimming outliers
#' @param showticks Whether to show axis ticks
#' @param shownumbers Whether to show axis numbers
#' @param showlabels Whether to show axis labels
#' @param plot.margins plotting margins
#' @param res Resolution for plotting
#' @param dens_res Resolution for density estimation
#' @param population_colors Colors for populations
#' @param auto_assign_color Whether to auto-assign colors
#' @param background_color Color for background
#' @param rgba_blending Whether to use RGBA blending
#' @param interpolate Whether to interpolate plot
#' @param rare_on_top Whether to plot rare populations on top
#' @param highlight_pops Populations to highlight
#' @param transformation Transformation type (e.g., arcsinh, logicle)
#' @param tr.parameters Parameters for transformations
#' @param print_estimated_parameters Whether to print estimated parameters
#' @param default_b Default value for arcsinh transformation
#' @param log_fsc Whether to log-transform FSC channels
#' @param log_ssc Whether to log-transform SSC channels
#' @param label_size Size for axis labels
#' @param color_pops Whether to enable population coloring
#' @param gates Gates to be drawn
#' @param old_method Whether to use older rendering method
#' @param density.overlay Whether to overlay density plots
#' @param max_points Maximum number of points to plot (for downsampling)
#' @param fulllength_label binary, whether to give full length label,default T
#'
#' @return A plot of the flow cytometry data
#' @export
PlotFlowFrame <- function(
    obj, # FlowFrame, GatingSet, or FlowSet object
    channels, # Channels to be plotted
    colormap_marker = NA,
    pop_idex.list = NA, # Index list of populations
    plotting_background = FALSE, # Plot background data
    title, # Plot title
    title.font.size, # Font size for title
    subtitle,
    subtitle.color,
    comp = FALSE, # Apply compensation (TRUE/FALSE)
    spillovermatix = NA, # Spillover matrix for compensation
    min.xy, # Minimum values for x and y axes
    max.xy,
    lowerlimit.quantile = 0.0003, # Lower limit quantile for trimming outliers
    showticks = TRUE, # Show axis ticks
    shownumbers = FALSE, # Show axis numbers
    showlabels = TRUE, # Show axis labels
    plot.margins = c(2, 2, 2, 2), # margin
    res = NA, # Resolution for plotting
    dens_res, # Resolution for density estimation
    population_colors = Color103, # Colors for populations
    auto_assign_color = TRUE, # Auto-assign color
    background_color = "#525252",
    margin_color,
    rgba_blending = TRUE, # Apply RGBA blending
    interpolate = FALSE, # Interpolate plot
    rare_on_top = TRUE, # Plot rare populations on top
    highlight_pops = NA, # Highlight specific populations
    highlight.alpha = 0.6,
    transformation = "logicle", # Transformation type (e.g., arcsinh, logicle)
    tr.parameters = NA, # Parameters for transformations
    print_estimated_parameters = FALSE, # Print estimated transformation parameters
    default_b = 400, # Default value for arcsinh transformation
    log_fsc = FALSE, # Log-transform FSC channels
    log_ssc = FALSE, # Log-transform SSC channels
    label_size = 0.8, # Axis label size
    color_pops = TRUE, # Enable population coloring
    gates, # Gates to be drawn
    gates_color = "grey", # Gates color
    old_method = TRUE, # Use older rendering method
    density.overlay = FALSE,
    max_points = 500000, # Maximum number of points to plot
    fulllength_label = T
    ) {
  # Section 1: Object Handling -------------------------------
  # Ensure `obj` is a valid FlowFrame or extract the florframe from root node from GatingSet
  if (class(obj) == "flowFrame") {
  } else if (class(obj) == "cytoframe") {
    obj <- cytoframe_to_flowFrame(obj)
  } else { # handling sets
    message("GatingSet or flowSet, only plot the 1st element")
    if (class(obj) == "GatingSet") {
      obj <- gs_pop_get_data(obj, "root") # get the root
    }
    if (class(obj) == "flowSet") {
      obj <- obj
    }

    obj <- obj[[1]] # a flowset
  }
  if(missing(title)) {title = identifier(obj)}

  # Section 2: Marker and Channel Validation -----------------
  # Parse channels and ensure they are valid in `obj`
  channels <- f.parseChannels(obj, channels)


  if (length(channels) < 2) {
    stop("At least one marker is not presented in the data!")
  }
  #print(parameters(obj)[["name"]])
  p_chnames <- parameters(obj)[["name"]][channels] # somehow using colnames(obj) is not working
  p_markers <- as.vector(parameters(obj)[["desc"]])[channels]
  #print(colnames(obj))
  #print(p_chnames)
  # Ensure p_chnames and p_markers are not NULL or empty
  if (length(p_chnames) == 0 || is.null(p_chnames)) {
    stop("Channel names could not be determined. Check that the channels parameter is valid.")
  }

  if (length(p_markers) == 0) {
    p_markers <- rep("", length(p_chnames))
  }

  # test if its plotting fsc/ssc channel
  if (sum(grepl("SC|ime", p_markers)) == 2) {
    comp <- FALSE
  } # switch off comp if both are linear
  p_markers[is.na(p_markers)] <- "" # for those markers that dont have desc

  # Section 3: Compensation Logic ----------------------------
  # Apply compensation matrix if requested
  if (comp == TRUE) {
    if (sum(grepl("Comp", p_chnames)) == 0) { # detect if there is Comp in the tittle
      if (is.na(spillovermatix)) {
        spillovermatix <- spillover(obj)
        spillovermatix <- spillovermatix[!unlist(lapply(spillover(obj), is.null))][[1]]
      }
      obj <- compensate(obj, spillovermatix)
      label_prefix <- "Comp-"
    } else {
      p_chnames <- gsub("Comp-", "", p_chnames)
      label_prefix <- "Comp-"
    }
  } else {
    label_prefix <- ""
  }

  # Section 4: Population Index Handling ---------------------
  # Subset populations to plot based on provided indices
  if (missing(pop_idex.list) | length(na.omit(pop_idex.list)) == 0) {
    # pop_idex.list <-NA
  } else {
    #' update the population info
    #' if not population color is supplied, using color103
    population_colors_updated <- lapply(seq(length(pop_idex.list)), function(x) {
      if (!names(pop_idex.list)[x] %in% names(population_colors)) { # if the color is not supplied in the population_colors
        # autoassign color could be useful in unsupervised analysis, but in supervised analysis such as flowdensity based gating
        # it could be confusing
        if (auto_assign_color == TRUE) {
          warning(paste0("No color assined for ", names(pop_idex.list)[x], ", using  color ", Color103[x]))
          Color103[x]
        } else {
          background_color
        }
      } else {
        population_colors[[names(pop_idex.list)[x]]]
      }
    })

    names(population_colors_updated) <- names(pop_idex.list)
    population_colors_updated <- unlist(population_colors_updated)

    if (plotting_background == TRUE) { # add "root" to the pop_idex.list
      r <- setdiff(seq(nrow(obj)), unlist(pop_idex.list))
      if (length(r) > 0) {
        pop_idex.list <- c(pop_idex.list, list("root" = r)) # the root here has already removed those foreground data pointts
        population_colors_updated <- c(population_colors_updated, background_color)
        names(population_colors_updated) <- names(pop_idex.list)

      }
    }
  }
  #print(population_colors_updated)
  #print(pop_idex.list)

  # get the data
  dat <- exprs(obj)[, channels]

  # Downsample if needed
  if(nrow(dat) > max_points) {
    set.seed(42)  # For reproducibility
    idx <- sample(nrow(dat), max_points)
    dat <- dat[idx, ]

    # If using population indices, also subset those
    if(!missing(pop_idex.list) && length(na.omit(pop_idex.list)) > 0) {
      pop_idex.list <- lapply(pop_idex.list, function(x) {
        intersect(x, idx)
      })
    }
  }

  # Section 5: Data Transformation ---------------------------
  # Apply transformation (e.g., arcsinh, logicle) based on input
  if (!is.null(obj@parameters@data$trans)) {
    trans_parameters_used <- data.frame(
      row.names = obj@parameters@data$name,
      trans = obj@parameters@data$trans
    )
  } else {
    trans_parameters_used <- NULL
  }

  ranges <- list(c(), c())
  major_tick_loc <- list(c(), c())
  minor_tick_loc <- list(c(), c())

  if (missing(min.xy)) {
    min.xy <- c(NA, NA)
  } else {
    if (length(min.xy) == 1) {
      min.xy <- c(min.xy[1], min.xy[1])
    }
  }

  if (missing(max.xy)) {
    max.xy <- c(NA, NA)
  } else {
    if (length(max.xy) == 1) {
      max.xy <- c(max.xy[1], max.xy[1])
    }
  }

  for (i in 1:2) {
    # Ensure p_chnames[i] exists before using it in grepl
    if (i > length(p_chnames)) {
      stop(paste("Channel index", i, "is out of bounds. Only", length(p_chnames), "channels available."))
    }

    #' Transform the linear channels
    if (grepl("SC", p_chnames[i])) {
      # if its FSC or SSC, assin the min
      if (grepl("SC", p_chnames[i])) {
        min.xy[i] <- 30
      } else {
        min.xy[i] <- as.numeric(parameters(obj)[["minRange"]][channels[i]])
      }

      if (is.na(max.xy[i])) {
        plot_ceiling <- as.numeric(parameters(obj)[["maxRange"]][channels[i]])
      } else {
        plot_ceiling <- max.xy[i]
      }

      if (grepl("FSC", p_chnames[i])) {
        if (log_fsc == TRUE) {
          mtd <- "log"
        } else {
          mtd <- "linear"
        }
      }
      if (grepl("SSC", p_chnames[i])) {
        if (log_ssc == TRUE) {
          min.xy[i] <- 3000
          plot_ceiling <- 262144
          mtd <- "log"
        } else {
          mtd <- "linear"
        }
      }

      lim <- c(min.xy[i], plot_ceiling)
      dat[, i][dat[, i] < min.xy[i]] <- min.xy[i]
      dat[, i][dat[, i] > plot_ceiling] <- plot_ceiling

      ranges[[i]] <- f.transFscSscData(lim, method = mtd)
      dat[, i] <- f.transFscSscData(dat[, i], method = mtd)

      if (mtd == "linear") {
        major_tick_loc[[i]] <- major_tick_linear
        minor_tick_loc[[i]] <- minor_tick_linear
      } else {
        major_tick_loc[[i]] <- f.transFscSscData(major_tick, method = mtd)
        minor_tick_loc[[i]] <- f.transFscSscData(minor_tick, method = mtd)
      }
    } else if (grepl("ime", p_chnames[i])) {
      ranges[[i]] <- c(0, max(dat[, i]))
      time_ticks <- setNames(seq(0, max(dat[, i]), 1000), seq(0, max(dat[, i]), 1000))
      major_tick_loc[[i]] <- time_ticks
      minor_tick_loc[[i]] <- NA
    } else if (grepl("PC|UMAP", p_chnames[i])) {
      if (is.na(min.xy[i])) {
        min.xy[i] <- as.numeric(parameters(obj)[["minRange"]][channels[i]])
      }

      if (is.na(max.xy[i])) {
        plot_ceiling <- as.numeric(parameters(obj)[["maxRange"]][channels[i]])
      } else {
        plot_ceiling <- max.xy[i]
      }

      lim <- c(min.xy[i], plot_ceiling)
      dat[, i] <- pmax(dat[, i], lim[1])
      dat[, i] <- pmin(dat[, i], lim[2])

      ranges[[i]] <- f.transFscSscData(lim, method = "linear")
      dat[, i] <- f.transFscSscData(dat[, i], method = "linear")

      major_tick_loc[[i]] <- seq(-10, 10, 1)
      minor_tick_loc[[i]] <- seq(-10, 10, 0.5)
    } else {
      # Transformation the Fluorescence channels
      ## set the minimun plotting value to 30 for log trans to "hide" the "fences"
      if (transformation == "log") {
        min.xy[i] <- 30
      }

      # using maxRnage and channles info, so the plot ceiling could be more specific
      if(is.na(min.xy[i])) {
        plot_floor <- parameters(obj)[["minRange"]][channels[i]]
      } else {
        plot_floor <-min.xy[i]
      }

      if(is.na(max.xy[i])) {
        plot_ceiling <- parameters(obj)[["maxRange"]][channels[i]]
      } else {
        plot_ceiling <-max.xy[i]
      }

      #' if the maximun is less then 10, the data is probably transformed
      if (plot_ceiling < 10) {
        # condition that the data is already transformed
        if (is.na(min.xy[i]) | (min.xy[i] < -10)) {
          bel <- round(sum(dat[, i] <= plot_floor) / nrow(dat) * 100, 1)
          if (bel >= 10) {
            warning(paste0(bel, " % of events on axis"))
          }

          lim <- c(plot_floor, plot_ceiling)
          dat[, i][dat[, i] < lim[1]] <- lim[1]
        } else {
          lim <- c(min.xy[i], plot_ceiling)
          dat[, i][dat[, i] < lim[1]] <- lim[1]
        }

        ranges[[i]] <- lim
        if (is.null(trans_parameters_used)) {
          # Showing transformed value
          major_tick_loc[[i]] <- unique(c(ceiling(lim[1]), c(0, 1, 2, 3, 4)))
          names(major_tick_loc[[i]]) <- major_tick_loc[[i]] # cautious here
          minor_tick_loc[[i]] <- major_tick_loc[[i]] + 0.5
          names(minor_tick_loc[[i]]) <- minor_tick_loc[[i]]
        } else {
          tr.parameter_used <- setNames(trans_parameters_used[p_chnames[i], ], "w")
          major_tick_loc[[i]] <- f.transVectData(major_tick, method = "logicle", tr.parameter_used)
          names(major_tick_loc[[i]]) <- major_tick # cautious here

          minor_tick_loc[[i]] <- f.transVectData(minor_tick, method = "logicle", tr.parameter_used)
          names(minor_tick_loc[[i]]) <- minor_tick
        }
      } else { # print("fluoresence channel")
        # Condition that the data is not transformed, which is the usual case
        # if no transformation parameters were supplied
        if (length(na.omit(tr.parameters)) == 0) {
          if (transformation == "arcsinh") {
            tr.parameter <- setNames(400, "b")
            message(
              paste0(
                "No transformation parameters, perform arcsinh transformation with b=400 for ",
                p_chnames[i]
              )
            )
          } else if (transformation == "logicle") {
            lgcl <- estimateLogicle(obj, channels = p_chnames[i]) # just estmate one channel instead of 2 so you dont have to estimate two times
            lgcl_unlisted <- unlist(summary(lgcl))
           # print(lgcl_unlisted)

            if(class(lgcl_unlisted) == "list") {

              tr.parameter <- round(c(lgcl_unlisted[[1]], lgcl_unlisted[[3]], lgcl_unlisted[[4]], lgcl_unlisted[[6]]),2)

            }else {

              tr.parameter <- round(lgcl_unlisted, 2)
            }
            names(tr.parameter) <- c("a", "m", "t", "w")

            message(
              paste0(
                "No transformation parameters, perform logicle transformation with estimated w=", tr.parameter["w"], " for ",
                p_chnames[i]
              )
            )
          } else {}
        } else { # do supply tr.parameters
          if (length(tr.parameters) == 1) {
            tr.parameters <- c(tr.parameters, tr.parameters)
          }

          if (class(tr.parameters) == "list") {
            tr.parameter <- tr.parameters[[p_chnames[i]]]
          } else if (class(tr.parameters) == "numeric") {
            if (transformation == "arcsinh") {
              if (tr.parameters[i] < 100) {
                message("b is too small, set b = 400")
                tr.parameter <- setNames(400, "b")
              } else {
                tr.parameter <- setNames(tr.parameters[i], "b")
              }
            } else if (transformation == "logicle") {
              if (tr.parameters[i] > 5) {
                message("w is too big, set w=0.8")
                tr.parameter <- setNames(0.8, "w")
              } else {
                tr.parameter <- setNames(tr.parameters[i], "w")
              }
            }
          }
        }
        # set minimal
        if (is.na(min.xy[i])) {
          if (is.na(lowerlimit.quantile)) {
            plot_floor <- min(dat[, i])
          } else {
            plot_floor <- quantile(dat[, i], lowerlimit.quantile)
          }
          lim <- c(plot_floor, plot_ceiling)
          dat[, i] <- pmax(dat[, i], lim[1])
        } else {
          lim <- c(min.xy[i], plot_ceiling)
          dat[, i] <- pmax(dat[, i], lim[1])
        }
        ## Transformation
        if (print_estimated_parameters) {
          print(unlist(tr.parameter))
        }

        ranges[[i]] <- f.transVectData(lim, method = transformation, tr.parameter)
        dat[, i] <- f.transVectData(dat[, i], method = transformation, tr.parameter)

        major_tick_loc[[i]] <- f.transVectData(major_tick, method = transformation, tr.parameter)
        names(major_tick_loc[[i]]) <- major_tick # cautious here

        minor_tick_loc[[i]] <- f.transVectData(minor_tick, method = transformation, tr.parameter)
        names(minor_tick_loc[[i]]) <- minor_tick
      }
    }
  }
  names(ranges) <- p_chnames

  if (is.na(res)) {
    if (missing(pop_idex.list)) {
      cn <- nrow(dat)
    } else {
      cn <- length(unlist(pop_idex.list))
    }

    if (cn >= 50000) {
      res <- 256
    } else if (cn < 50000 & cn >= 10000) {
      res <- 196
    } else if (cn < 10000 & cn >= 2500) {
      res <- 128
    } else if (cn < 2500 & cn >= 200) {
      res <- 96
    } else if (cn < 200 & cn >= 50) {
      res <- 64
    } else {
      res <- 32
    }
  } else {
    cn <- nrow(dat)
    res <- res
  }

  if (missing(dens_res)) {
    dens_res <- res
  }

  # Section 6: Plot Rendering -------------------------------
  # Generate and render the plot with specified options
  if (is.null(color_pops)) {
    color_pops <- FALSE
  }

  if (cn != 0) {
    if (!is.na(colormap_marker)) { # Heatmap
      colormap_idx <- f.parseChannels(obj, colormap_marker)
      color_name <- as.vector(parameters(obj)[["desc"]])[colormap_idx]
      # Use vectorized operation instead of apply
      na <- rowSums(is.na(dat)) > 0
      dat <- dat[!na, ]

      mk <- exprs(obj)[!na, colormap_idx]
      col <- number_to_color(mk, palette = "inferno")
      col <- col2rgb(col, alpha = TRUE)

      rgbwt <- scatter_points_rgbwt(dat,
        out_size = c(res, res),
        RGBA = col,
        xlim = ranges[[1]],
        ylim = ranges[[2]]
      )

      rgbwt <- rgbwt_to_rgba_int(rgbwt)
    } else {
      if (length(na.omit(pop_idex.list)) == 0 | color_pops == FALSE) {
        # No population info or no population color infor provided
        # Single ccolor will be used
        if (length(na.omit(pop_idex.list)) != 0) {
          l.dat <- lapply(pop_idex.list, function(x) {
            dat[x, ]
          })
          dat <- do.call(rbind, l.dat)
        }

        # Use vectorized operation instead of apply
        na <- rowSums(is.na(dat)) > 0
        dat <- dat[!na, ]
        # Plotting single population

        if (old_method) {
          rgbwt <- scatter_points_rgbwt(dat,
            out_size = c(res, res),
            RGBA = col2rgb("#25252505", alpha = TRUE),
            xlim = ranges[[1]],
            ylim = ranges[[2]]
          )

          rgbwt[, , 5] <- log1p(1 - rgbwt[, , 5]) / log(2) * 0.95
          rgbwt[, , 5][rgbwt[, , 5] == 0] <- 1
        } else {
          rgbwt <- scatter_points_rgbwt(dat,
            out_size = c(res, res),
            xlim = ranges[[1]],
            ylim = ranges[[2]]
          )
          rgbwt <- scale_toward_white(rgbwt, base_color = c(0, 0, 0), whiteness = 0.95)
        }

        rgbwt <- rgbwt_to_rgba_int(rgbwt)
      } else {
        # Plotting multiple populations
        l.dat <- lapply(pop_idex.list, function(x) {
          dat[x, , drop = FALSE]
        })
        #print(names(l.dat))
        # Pre-calculate population sizes
        pop_sizes <- vapply(l.dat, nrow, numeric(1))
        #print(pop_sizes)
        # change the list order to plot
        if (rare_on_top == TRUE) {
          l.dat <- l.dat[names(sort(pop_sizes, decreasing = FALSE))]
        }

        # plotting multiple populations by merging layers
        if (old_method) {
          l.rgbwt <- lapply(names(l.dat), function(x) {
            dat_subset <- l.dat[[x]]
            color <- population_colors_updated[x]
            color <- paste0(color, "10")
            rgbwt <- scatter_points_rgbwt(dat_subset,
              out_size = c(res, res),
              RGBA = col2rgb(color, alpha = TRUE),
              xlim = ranges[[1]],
              ylim = ranges[[2]]
            )

            if (x %in% highlight_pops) {
              rgbwt <- apply_kernel_rgbwt(rgbwt, "square", radius = 1)
            }

            rgbwt[, , 5] <- log1p(1 - rgbwt[, , 5]) / 0.8
            rgbwt[, , 5][rgbwt[, , 5] == 0] <- 1

            if (rgba_blending == TRUE) {
              rgbwt <- rgbwt_to_rgba_float(rgbwt)
            }

            return(rgbwt)
          })
        } else {

          l.rgbwt <- lapply(names(l.dat), function(x) {
            dat_subset <- l.dat[[x]]
            color <- population_colors_updated[x]
            rgbwt <- scatter_points_rgbwt(dat_subset,
              out_size = c(res, res),
              xlim = ranges[[1]],
              ylim = ranges[[2]]
            )
            if (x %in% highlight_pops) {
              cn <- nrow(dat_subset)

              if (cn >= 5000) {
                rds <- 1
              } else if (cn < 5000 & cn >= 2500) {
                rds <- 1
              } else if (cn < 2500 & cn >= 200) {
                rds <- 2
              } else {
                rds <- 3
              }

              rgbwt <- apply_kernel_rgbwt(rgbwt, "square", radius = rds)
            }

            if(nrow(dat_subset) > 5) {
              rgbwt <- scale_toward_white(rgbwt, base_color = c(col2rgb(color, alpha = FALSE)), whiteness = 0.95)
            }

            if (rgba_blending == TRUE) {
              rgbwt <- rgbwt_to_rgba_float(rgbwt)
            }

            return(rgbwt)
          })
          #print(length(l.rgbwt))
        }
        if (rgba_blending == TRUE) {
          if (old_method) {
            rgbwt <- blend_rgba_float(l.rgbwt)
          } else {

            hl_pop <- names(l.dat)[names(l.dat) %in% highlight_pops]

            if(length(hl_pop) >0 ) {
              hl_index <- which(names(l.dat) %in% hl_pop)
              l.rgbwt <- c(l.rgbwt[hl_index], l.rgbwt[-hl_index])


              rgbwt <- blend_rgba_float3(l.rgbwt, alpha = highlight.alpha) # no transparency if there is highlight pop

            }else {
              rgbwt <- blend_rgba_float3(l.rgbwt, alpha = 0.6) # it doesnt affect non overlapping pixels

            }



          }
          rgbwt <- rgba_float_to_rgba_int(rgbwt)
        } else {
          rgbwt <- merge_rgbwt(l.rgbwt)
          rgbwt <- rgbwt_to_rgba_int(rgbwt)
        }
      }
    }

    rstr <- rgba_int_to_raster(rgbwt)
  }

  # Section 7: Finalize plot -------------------------------
  # Handling plot ranges, ticks...
  usr <- c(ranges[[1]], ranges[[2]])
  par(

    pty = "s", bty = "o",
    mgp = c(3, 0.5, 0),
    mar = plot.margins#c(2, 2, 3, 2)
  )

  if(!missing(margin_color)){
    par(bg = margin_color)
  }

  plot(c(),
    xlim = ranges[[1]], xlab = "", xaxt = "n",
    ylim = ranges[[2]], ylab = "", yaxt = "n",
  )

  if(missing(title.font.size)) {
    if(nchar(title) > 30) {
      title.font.size = 0.6
    } else {
      title.font.size = 1
    }
  }

  title(title, adj = 0, line = 0.3, cex.main = title.font.size)



  if(missing(subtitle)) {
    if (!is.na(colormap_marker)) {
      title(color_name, adj = 1, line = 0.3, cex.main = title.font.size)
    }
  } else {
    if(missing(subtitle.color)){
      title(subtitle, adj = 1, line = 0.3, cex.main = title.font.size)
    } else {
      title(subtitle, adj = 1, line = 0.3, cex.main = title.font.size, col.main = subtitle.color)

    }

  }


  if(!missing(margin_color)){
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col="white", border= "black")
  par(new=TRUE)
  }

  if (cn != 0) {
    rasterImage(rstr,
      xleft = usr[1],
      xright = usr[2],
      ybottom = usr[3],
      ytop = usr[4],
      interpolate = interpolate
    )
  }

  # plot ticks and numbers
  if (showticks == TRUE) {
    f.axis_ticks(
      p_chnames, major_tick_loc, minor_tick_loc, ranges,
      shownumbers
    )
  }

  # add labels
  if (showlabels == TRUE) {
    if (showticks == TRUE & shownumbers == TRUE) {
      label_padding <- 1.3
    } else if (showticks == TRUE & shownumbers == FALSE) {
      label_padding <- 0.9
    } else {
      label_padding <- 0.5
    }
    axis_labels <- p_chnames

    if(fulllength_label == T) {
      # labels with both channel names and desciptions
      axis_labels[!grepl("SC|ime", p_chnames)] <-
        paste0(
          label_prefix, p_chnames[!grepl("SC|ime", p_chnames)],
          "::", p_markers[!grepl("SC|ime", p_chnames)]
        )
    } else {
      # labels with only descriptions
      axis_labels[!grepl("SC|ime", p_chnames)] <-
        paste0(
          label_prefix,  p_markers[!grepl("SC|ime", p_chnames)]
        )
    }

    if(log_ssc == T) {
      axis_labels[grepl("SSC", p_chnames)] <-
      paste0(
        p_chnames[grepl("SSC", p_chnames)], " LOG"
      )
    }

    if(log_fsc == T) {
      axis_labels[grepl("FSC", p_chnames)] <-
      paste0(
        p_chnames[grepl("FSC", p_chnames)], " LOG"
      )
    }

    f.axis_label(axis_labels, padding = label_padding, font.size = label_size)
  }




  if (!missing(gates) #| !is.na(gates)
  ) {
    f.plotGate(gates, xlim = ranges[[1]], ylim = ranges[[2]], log_ssc = log_ssc, col = gates_color)
  }

  # add gatings if neceessary


  if (density.overlay & cn > 10) {
    dens.color <- "red"
    dens.alpha = 0.5

    adjust.dens = 1
    x.dens <- density(dat[, 1], adjust = adjust.dens)
    y.dens <- density(dat[, 2], adjust = adjust.dens)

    x.axis <- x.dens$x#rescale(x.dens$x, xlim)
    y.axis <- y.dens$x#rescale(y.dens$y, ylim)
    x.dens <- rescale(x.dens$y, ranges[[2]])
    y.dens <- rescale(y.dens$y, ranges[[1]])


    lines(x.axis, x.dens, col = alpha(dens.color, dens.alpha) )
    lines(y.dens, y.axis, col = alpha(dens.color, dens.alpha))


  }



















}
