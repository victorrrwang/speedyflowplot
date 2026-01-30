# speedyflowplot 1.15

## New Features

* `flowPlot()` now supports multiple populations via the `ids` parameter
  * Multiple populations can be visualized with different colors
  * Uses alpha blending to combine overlapping populations
  * Color parameter accepts a vector of colors for each population

# speedyflowplot 1.14

## Package Structure Updates

* Updated package documentation to credit the scattermore package
* Added acknowledgments section in README files
* Updated DESCRIPTION to include scattermore attribution
* Added AI-generated content disclaimer to main README

# speedyflowplot 1.1

## New Features

* Added optimized functions using raw uint8_t arrays for improved performance:
  * `scatter_points_optimized()`: Optimized version of `scatter_points_rgbwt()`
  * `apply_kernel_optimized()`: Optimized version of `apply_kernel_rgbwt()`
  * `raw_to_raster()`: Convert raw uint8_t arrays to raster objects
  * `PlotFlowFrame_optimized()`: Optimized version of `PlotFlowFrame()`

* Added placeholder GPU functions (currently falling back to optimized CPU implementation):
  * `gpu_is_available()`: Check if GPU is available for acceleration
  * `scatter_points_gpu()`: GPU-accelerated version of `scatter_points_optimized()`
  * `apply_kernel_gpu()`: GPU-accelerated version of `apply_kernel_optimized()`
  * `PlotFlowFrame_gpu()`: GPU-accelerated version of `PlotFlowFrame_optimized()`

* Added documentation for optimization approach in `README_optimization.md`

* Added example scripts:
  * `examples/optimized_example.R`: Demonstrates performance improvements with optimized functions
  * `examples/gpu_example.R`: Demonstrates placeholder GPU functions

* Added documentation for potential GPU acceleration in `docs/gpu_acceleration.md`

## Performance Improvements

* Improved performance by using the `.Call()` interface instead of `.C()`
* Reduced memory usage by avoiding unnecessary data copies
* Optimized data transfer between R and C
* Direct manipulation of raw arrays in C for faster processing

## Technical Details

* Added new C++ implementation in `src/optimized_scatter.cpp`
* Updated `src/speedyflowplot.c` to register new functions
* Added R wrapper functions in `R/optimized_scatter.R`
* Added optimized main function in `R/PlotFlowFrame_optimized.R`
