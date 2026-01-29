# speedyflowplot

A fast plotting package for flow cytometry data visualization.

## Overview

speedyflowplot is a lightweight R package that provides efficient rendering of large scatterplots. It is particularly useful for flow cytometry data visualization where millions of points need to be plotted quickly.

This package contains core functions for flow cytometry visualization:

### Core Rendering Functions
- `scatter_points_rgbwt`: Render colored points into a RGBWT bitmap
- `rgbwt_to_rgba_int`: Convert a RGBWT matrix to an integer RGBA matrix
- `rgba_int_to_raster`: Create a raster from the given RGBA matrix
- `apply_kernel_rgbwt`: Apply a kernel to the given RGBWT raster
- `scale_toward_white`: Scale RGB colors toward white based on density
- `rgbwt_to_rgba_float`: Convert RGBWT matrix to floating-point RGBA matrix
- `rgba_float_to_rgba_int`: Convert a float RGBA bitmap to integer RGBA bitmap
- `blend_rgba_float`: Blend RGBA matrices
- `blend_rgba_float3`: Blend multiple RGBA arrays with selective transparency
- `merge_rgbwt`: Merge RGBWT matrices

### Flow Cytometry Plotting Functions
- `flowPlot`: Simple 2D density plot for flow cytometry data
- `PlotFlowFrame`: Advanced flow cytometry plotting with multiple options

### Utility Functions
- `f.parseChannels`: Parse channel names or indices from a flow cytometry object
- `f.transVectData`: Transform flow cytometry data
- `f.transFscSscData`: Transform FSC/SSC data
- `f.axis_ticks`: Draw axis ticks for flow cytometry plots
- `f.axis_label`: Add axis labels to flow cytometry plots
- `f.plotGate`: Plot gates on flow cytometry plots
- `number_to_color`: Map numeric values to color gradients

## Installation

```r
# Install from GitHub
devtools::install_github("yourusername/speedyflowplot")
```

## Usage

### Basic Usage

```r
library(speedyflowplot)

# Create a simple scatterplot
xy <- cbind(rnorm(1e6), rnorm(1e6))
rgba <- col2rgb(rgb(0.25, 0.5, 0.75, 0.04), alpha = TRUE)
rgbwt <- scatter_points_rgbwt(xy, rgba = rgba)
rgba_int <- rgbwt_to_rgba_int(rgbwt)
raster <- rgba_int_to_raster(rgba_int)

# Plot the raster
plot(c(), xlim = c(-3, 3), ylim = c(-3, 3))
rasterImage(raster, -3, -3, 3, 3)
```

### Flow Cytometry Plotting

```r
# Simple flow cytometry plot
flowPlot(flowframe_obj, channels = c("FSC-A", "SSC-A"), 
         title = "Flow Cytometry Plot", 
         density.overlay = TRUE)

# Advanced flow cytometry plot
PlotFlowFrame(flowframe_obj, channels = c("FSC-A", "SSC-A"),
              transformation = "logicle",
              showticks = TRUE, shownumbers = TRUE)
```

## Acknowledgments

This package incorporates code and concepts from the [**scattermore**](https://github.com/exaexa/scattermore) package by Mirek Kratochvil and Tereza Kulichova. The core C++ rendering functions for scatter plot rasterization are based on scattermore's implementation.

## License

GPL (>= 3)
