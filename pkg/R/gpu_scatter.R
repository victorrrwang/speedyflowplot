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

#' Check if GPU is available
#'
#' This function checks if a compatible GPU is available for acceleration.
#' Currently, this is a placeholder function that always returns FALSE.
#' In a full implementation, this would check for OpenCL or CUDA compatibility.
#'
#' @return A logical value indicating whether a compatible GPU is available.
#'
#' @export
gpu_is_available <- function() {
  # This is a placeholder function
  # In a real implementation, this would check for OpenCL or CUDA compatibility
  message("GPU acceleration is not yet implemented. This is a placeholder function.")
  return(FALSE)
}

#' Scatter points using GPU acceleration
#'
#' This function is a placeholder for GPU-accelerated scatter points rendering.
#' Currently, it falls back to the optimized CPU implementation.
#' In a full implementation, this would use OpenCL or CUDA for GPU acceleration.
#'
#' @inheritParams scatter_points_optimized
#' @param force_gpu Logical value indicating whether to force GPU usage even if not recommended.
#'
#' @return A raw vector containing RGBA bitmap data.
#'
#' @export
scatter_points_gpu <- function(xy,
                              xlim = c(min(xy[, 1]), max(xy[, 1])),
                              ylim = c(min(xy[, 2]), max(xy[, 2])),
                              out_size = c(512, 512),
                              RGBA = c(0, 0, 0, 255),
                              map = NULL,
                              palette = NULL,
                              force_gpu = FALSE) {
  # Check if GPU is available
  if (!gpu_is_available() && !force_gpu) {
    message("GPU not available or not recommended. Falling back to optimized CPU implementation.")
    return(scatter_points_optimized(xy, xlim, ylim, out_size, RGBA, map, palette))
  }
  
  # This is a placeholder for GPU implementation
  # In a real implementation, this would use OpenCL or CUDA with .Call to C++ functions
  # For example:
  # raw_output <- raw(out_size[1] * out_size[2] * 4)
  # .Call("scatter_points_gpu", xy, xlim, ylim, as.integer(out_size), RGBA, raw_output, PACKAGE = "speedyflowplot")
  
  message("GPU acceleration is not yet implemented. Falling back to optimized CPU implementation.")
  return(scatter_points_optimized(xy, xlim, ylim, out_size, RGBA, map, palette))
}

#' Apply kernel using GPU acceleration
#'
#' This function is a placeholder for GPU-accelerated kernel application.
#' Currently, it falls back to the optimized CPU implementation.
#' In a full implementation, this would use OpenCL or CUDA for GPU acceleration.
#'
#' @inheritParams apply_kernel_optimized
#' @param force_gpu Logical value indicating whether to force GPU usage even if not recommended.
#'
#' @return A raw vector containing the blurred RGBA bitmap data.
#'
#' @export
apply_kernel_gpu <- function(raw_input,
                            out_size,
                            filter = "circle",
                            mask = default_kernel(filter, radius, sigma),
                            radius = 2,
                            sigma = radius / 2,
                            force_gpu = FALSE) {
  # Check if GPU is available
  if (!gpu_is_available() && !force_gpu) {
    message("GPU not available or not recommended. Falling back to optimized CPU implementation.")
    return(apply_kernel_optimized(raw_input, out_size, filter, mask, radius, sigma))
  }
  
  # This is a placeholder for GPU implementation
  # In a real implementation, this would use OpenCL or CUDA
  message("GPU acceleration is not yet implemented. Falling back to optimized CPU implementation.")
  return(apply_kernel_optimized(raw_input, out_size, filter, mask, radius, sigma))
}

#' Plot flow frame using GPU acceleration
#'
#' This function is a placeholder for GPU-accelerated flow frame plotting.
#' Currently, it falls back to the optimized CPU implementation.
#' In a full implementation, this would use OpenCL or CUDA for GPU acceleration.
#'
#' @inheritParams PlotFlowFrame_optimized
#' @param use_gpu Logical value indicating whether to use GPU acceleration.
#' @param force_gpu Logical value indicating whether to force GPU usage even if not recommended.
#'
#' @return A plot of the flow cytometry data
#'
#' @export
PlotFlowFrame_gpu <- function(
    obj,
    channels,
    colormap_marker = NA,
    pop_idex.list = NA,
    plotting_background = FALSE,
    title,
    title.font.size,
    comp = FALSE,
    spillovermatix = NA,
    min.xy,
    max.xy,
    lowerlimit.quantile = 0.0003,
    showticks = TRUE,
    shownumbers = FALSE,
    showlabels = TRUE,
    res = NA,
    dens_res,
    population_colors = autoflow_colors,
    auto_assign_color = TRUE,
    background_color = "#525252",
    rgba_blending = TRUE,
    interpolate = FALSE,
    rare_on_top = TRUE,
    highlight_pops = NA,
    transformation = "logicle",
    tr.parameters = NA,
    print_estimated_parameters = FALSE,
    default_b = 400,
    log_fsc = FALSE,
    log_ssc = FALSE,
    label_size = 0.8,
    gate_by_Density = FALSE,
    color_pops = TRUE,
    gates,
    old_method = TRUE,
    density.overlay = FALSE,
    max_points = 500000,
    use_raw_arrays = TRUE,
    use_gpu = TRUE,
    force_gpu = FALSE
    ) {
  
  # Check if GPU is available
  if (!gpu_is_available() && !force_gpu) {
    message("GPU not available or not recommended. Falling back to optimized CPU implementation.")
    return(PlotFlowFrame_optimized(
      obj = obj,
      channels = channels,
      colormap_marker = colormap_marker,
      pop_idex.list = pop_idex.list,
      plotting_background = plotting_background,
      title = title,
      title.font.size = title.font.size,
      comp = comp,
      spillovermatix = spillovermatix,
      min.xy = min.xy,
      max.xy = max.xy,
      lowerlimit.quantile = lowerlimit.quantile,
      showticks = showticks,
      shownumbers = shownumbers,
      showlabels = showlabels,
      res = res,
      dens_res = dens_res,
      population_colors = population_colors,
      auto_assign_color = auto_assign_color,
      background_color = background_color,
      rgba_blending = rgba_blending,
      interpolate = interpolate,
      rare_on_top = rare_on_top,
      highlight_pops = highlight_pops,
      transformation = transformation,
      tr.parameters = tr.parameters,
      print_estimated_parameters = print_estimated_parameters,
      default_b = default_b,
      log_fsc = log_fsc,
      log_ssc = log_ssc,
      label_size = label_size,
      gate_by_Density = gate_by_Density,
      color_pops = color_pops,
      gates = gates,
      old_method = old_method,
      density.overlay = density.overlay,
      max_points = max_points,
      use_raw_arrays = use_raw_arrays
    ))
  }
  
  # This is a placeholder for GPU implementation
  # In a real implementation, this would use OpenCL or CUDA
  message("GPU acceleration is not yet implemented. Falling back to optimized CPU implementation.")
  return(PlotFlowFrame_optimized(
    obj = obj,
    channels = channels,
    colormap_marker = colormap_marker,
    pop_idex.list = pop_idex.list,
    plotting_background = plotting_background,
    title = title,
    title.font.size = title.font.size,
    comp = comp,
    spillovermatix = spillovermatix,
    min.xy = min.xy,
    max.xy = max.xy,
    lowerlimit.quantile = lowerlimit.quantile,
    showticks = showticks,
    shownumbers = shownumbers,
    showlabels = showlabels,
    res = res,
    dens_res = dens_res,
    population_colors = population_colors,
    auto_assign_color = auto_assign_color,
    background_color = background_color,
    rgba_blending = rgba_blending,
    interpolate = interpolate,
    rare_on_top = rare_on_top,
    highlight_pops = highlight_pops,
    transformation = transformation,
    tr.parameters = tr.parameters,
    print_estimated_parameters = print_estimated_parameters,
    default_b = default_b,
    log_fsc = log_fsc,
    log_ssc = log_ssc,
    label_size = label_size,
    gate_by_Density = gate_by_Density,
    color_pops = color_pops,
    gates = gates,
    old_method = old_method,
    density.overlay = density.overlay,
    max_points = max_points,
    use_raw_arrays = use_raw_arrays
  ))
}
