# speedyflowplot

**Fast Flow Cytometry Plotting with C++ Acceleration**

This repository contains the **speedyflowplot** R package, a high-performance visualization tool for flow cytometry data. The package provides efficient rendering of large scatterplots (millions of points) through optimized C++ implementations with multi-threaded processing.


## Overview

speedyflowplot is an R package that combines R's ease of use with C++ performance for flow cytometry data visualization. The package implements a sophisticated rendering pipeline that converts scatter plot data into high-quality raster images with support for density visualization, color blending, and kernel-based smoothing.

### Key Features

- **High Performance**: C++ backend with optimized memory management and data transfer
- **Large Dataset Support**: Efficiently handles millions of data points
- **Multiple Rendering Modes**: Single color, multi-color, and indexed color rendering
- **Density Visualization**: Built-in support for density-based color scaling
- **Kernel Operations**: Blur and smoothing operations for enhanced visualization
- **Flow Cytometry Integration**: Native support for `flowCore` objects

## Architecture

### Core Rendering Pipeline

The package implements a multi-stage rendering pipeline:

```
Scatter Points → RGBWT Bitmap → RGBA Matrix → Raster Image
```

1. **Scatter Point Rendering**: Converts (x, y) coordinates into a RGBWT (Red, Green, Blue, Weight, Total) bitmap
2. **Color Conversion**: Transforms RGBWT format to RGBA (Red, Green, Blue, Alpha) matrices
3. **Raster Generation**: Creates R graphics raster objects for plotting

### Implementation Details

#### C++ Backend

The package uses two interfaces for C++ integration:

- **`.C()` Interface**: Original implementation for RGBWT-based rendering
  - `scatter_indexed_rgbwt`: Indexed color rendering
  - `scatter_singlecolor_rgbwt`: Single color rendering
  - `scatter_multicolor_rgbwt`: Multi-color rendering
  - `kernel_rgbwt`: Kernel-based blur operations

- **`.Call()` Interface**: Optimized implementation using raw uint8_t arrays
  - `scatter_points_raw`: Direct scatter point rendering to raw RGBA arrays
  - `scatter_indexed_raw`: Indexed rendering with raw arrays
  - `apply_kernel_raw`: Kernel operations on raw RGBA arrays

#### Performance Optimizations

The optimized functions (`*_optimized`) provide significant performance improvements through:

1. **Direct Memory Access**: Using `.Call()` interface to avoid data copying between R and C++
2. **Raw Array Operations**: Working directly with `uint8_t` arrays instead of converting between formats
3. **Reduced Conversions**: Eliminating intermediate RGBWT → RGBA → Raster conversions
4. **Multi-threading**: Parallel processing for large datasets (via thread blocks)

#### Memory Management

- **Zero-copy Operations**: Direct pointer access to R's internal data structures
- **Efficient Data Transfer**: Raw arrays minimize memory overhead
- **In-place Operations**: Kernel operations can modify arrays in-place when possible

## Package Structure

```
speedyflowplot/
├── pkg/                          # Main R package directory
│   ├── R/                        # R wrapper functions
│   │   ├── scatter_points_rgbwt.R      # Original RGBWT functions
│   │   ├── optimized_scatter.R         # Optimized raw array functions
│   │   ├── PlotFlowFrame.R             # Main plotting function
│   │   ├── PlotFlowFrame_optimized.R   # Optimized plotting function
│   │   └── ...
│   ├── src/                      # C++ source code
│   │   ├── optimized_scatter.cpp       # Optimized .Call() implementations
│   │   ├── scatter_*_rgbwt.cpp         # Original .C() implementations
│   │   ├── kernel_rgbwt.cpp            # Kernel operations
│   │   ├── scatters.h                  # Scatter function declarations
│   │   └── kernels.h                   # Kernel function declarations
│   ├── examples/                 # Example scripts
│   ├── man/                      # R documentation
│   └── DESCRIPTION               # Package metadata
└── releases/                     # Archived package versions
```

## Installation

### From Source

To install the package from source:

```r
# Install dependencies first
install.packages(c("grDevices", "graphics", "stats", "scales", "viridis"))

# If you don't have BiocManager installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install(c("flowCore", "flowWorkspace"))

# Install devtools if you don't have it
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install the package from the local directory
devtools::install("pkg")
```

Replace `"pkg"` with the actual path to the package directory if you're not in the repository root.

### Quick Installation

```r
# Install dependencies
install.packages(c("grDevices", "graphics", "stats", "scales", "viridis"))

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("flowCore")

# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install from local directory
devtools::install("pkg")
```

### Dependencies

**Required:**
- R (>= 3.5.0)
- `flowCore` (Bioconductor)
- `grDevices`, `graphics`, `stats` (base R)
- `scales`, `viridis` (CRAN)

**Optional:**
- `flowWorkspace` (Bioconductor) - for advanced flow cytometry operations

## Usage

### Basic Flow Cytometry Plotting

```r
library(speedyflowplot)
library(flowCore)

# Load flow cytometry data
ff <- read.FCS("your_data.fcs")

# Simple density plot
flowPlot(ff, 
         channels = c("FSC-A", "SSC-A"), 
         title = "Flow Cytometry Plot", 
         density.overlay = TRUE)

# Advanced plotting with transformations
PlotFlowFrame(ff, 
              channels = c("FSC-A", "SSC-A"),
              transformation = "logicle",
              showticks = TRUE, 
              shownumbers = TRUE)
```

### Optimized Functions

For better performance with large datasets, use the optimized functions. The optimized functions have the same parameters as the original functions, making them easy to use as drop-in replacements:

```r
# Optimized scatter point rendering
xy <- cbind(rnorm(1e6), rnorm(1e6))
raw_output <- scatter_points_optimized(
    xy = xy,
    out_size = c(512, 512),
    xlim = c(-3, 3),
    ylim = c(-3, 3)
)
raster <- raw_to_raster(raw_output, c(512, 512))

# Optimized flow cytometry plotting
PlotFlowFrame_optimized(
    obj = ff,
    channels = c("FSC-A", "SSC-A"),
    title = "My Plot"
)
```

### Using the GPU Functions (Placeholder)

The package includes placeholder GPU functions that currently fall back to the optimized CPU implementation:

```r
# Check if GPU is available (always returns FALSE in this version)
gpu_available <- gpu_is_available()

# Use the GPU version of PlotFlowFrame
PlotFlowFrame_gpu(
    obj = ff,
    channels = c("FSC-A", "SSC-A"),
    title = "GPU Plot"
)
```

### Performance Comparison

To compare the performance of the original and optimized functions:

```r
library(microbenchmark)

# Compare original vs optimized
benchmark <- microbenchmark(
    original = {
        rgbwt <- scatter_points_rgbwt(xy, out_size = c(512, 512))
        rgba <- rgbwt_to_rgba_int(rgbwt)
        raster <- rgba_int_to_raster(rgba)
    },
    optimized = {
        raw_output <- scatter_points_optimized(xy, out_size = c(512, 512))
        raster <- raw_to_raster(raw_output, c(512, 512))
    },
    times = 10
)

print(benchmark)

# Compare PlotFlowFrame functions
benchmark_plot <- microbenchmark(
    original = {
        PlotFlowFrame(
            obj = ff,
            channels = c("FSC-A", "SSC-A")
        )
    },
    optimized = {
        PlotFlowFrame_optimized(
            obj = ff,
            channels = c("FSC-A", "SSC-A")
        )
    },
    times = 5
)

print(benchmark_plot)
```

### Example Scripts

The package includes example scripts that demonstrate the new functions:

- `pkg/examples/optimized_example.R`: Demonstrates the optimized functions
- `pkg/examples/gpu_example.R`: Demonstrates the placeholder GPU functions

To run these examples:

```r
# Source the example scripts
source(system.file("examples", "optimized_example.R", package = "speedyflowplot"))
source(system.file("examples", "gpu_example.R", package = "speedyflowplot"))
```

See [`pkg/examples/optimized_example.R`](pkg/examples/optimized_example.R) for a complete performance comparison example.

## Technical Details

### Optimization Approach

The optimized functions use the `.Call()` interface instead of `.C()` to achieve better performance:

1. **No Data Copying**: Direct access to R's internal data structures via `SEXP` objects
2. **Raw Array Operations**: Working with `uint8_t` arrays eliminates format conversion overhead
3. **Efficient Memory Layout**: Contiguous memory arrays for better cache locality
4. **Reduced Function Calls**: Fewer intermediate R function calls in the rendering pipeline

For detailed information, see [`pkg/README_optimization.md`](pkg/README_optimization.md).

### C++ Implementation

The C++ code is organized into several modules:

- **Scatter Rendering**: Functions for converting point coordinates to bitmaps
- **Kernel Operations**: Convolution operations for blur and smoothing
- **Color Conversion**: Efficient conversion between color formats
- **Thread Management**: Multi-threaded processing support

Key C++ functions:

- `scatter_points_raw()`: Main optimized scatter rendering function
- `scatter_indexed_raw()`: Indexed color rendering with palette mapping
- `apply_kernel_raw()`: Kernel convolution on RGBA arrays

### Data Formats

- **RGBWT**: Red, Green, Blue, Weight, Total - intermediate format for density accumulation
- **RGBA**: Red, Green, Blue, Alpha - final color format with transparency
- **Raw Arrays**: `uint8_t` arrays for direct memory manipulation

## Documentation

- **Build Instructions**: [`build.md`](build.md) - How to build the package from source
- **Optimization Details**: [`pkg/README_optimization.md`](pkg/README_optimization.md) - Technical details on performance optimizations
- **Package README**: [`pkg/README.md`](pkg/README.md) - Package-specific documentation
- **R Documentation**: See `?speedyflowplot` after installation

## Version

Current version: **1.14**

See [`pkg/NEWS.md`](pkg/NEWS.md) for version history and changelog.

## Contributors

- **Victor Wang** - Author (victorrwang@gmail.com)

## Acknowledgments

This package incorporates code and concepts from the [**scattermore**](https://github.com/exaexa/scattermore) package by Mirek Kratochvil and Tereza Kulichova. The core C++ rendering functions for scatter plot rasterization are based on scattermore's implementation, which provides fast scatterplot rendering for R. We are grateful to the scattermore developers for their excellent work on high-performance scatterplot visualization.

**scattermore** is available on CRAN and GitHub:
- CRAN: `install.packages('scattermore')`
- GitHub: https://github.com/exaexa/scattermore
- Website: https://exaexa.github.io/scattermore/

## License

GPL (>= 3)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Note

This README file was generated with the assistance of AI. While efforts have been made to ensure accuracy, please verify any technical details and refer to the package source code and documentation for authoritative information.


