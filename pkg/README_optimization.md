# Optimizing speedyflowplot with Raw uint8_t Arrays

This document explains the optimization approach implemented in the speedyflowplot package to improve performance by using raw uint8_t arrays for data transfer between R and C.

## Background

The original implementation of speedyflowplot uses the `.C()` interface to call C functions from R, which is based on the [scattermore](https://github.com/exaexa/scattermore) package's implementation. This interface has several limitations:

1. It makes copies of R data when passing to C
2. It converts between R and C data types (e.g., R numeric to C float)
3. It requires copying data back to R

These limitations can lead to performance bottlenecks, especially when dealing with large datasets or high-resolution plots.

## Optimization Approach

The optimization approach implemented in this package uses the `.Call()` interface instead of `.C()`. The `.Call()` interface allows direct access to R's internal data structures without copying, which can significantly improve performance.

The key optimizations include:

1. **Using raw uint8_t arrays**: Instead of converting between different formats (RGBWT, RGBA, raster), we work directly with raw uint8_t arrays.
2. **Avoiding unnecessary data copies**: We pass pointers to data between R and C, avoiding unnecessary copies.
3. **Direct manipulation of raw arrays**: We manipulate the raw arrays directly in C, avoiding the overhead of R function calls.

## Implementation

The optimization is implemented in the following files:

- `src/optimized_scatter.cpp`: C++ implementation of the optimized functions
- `R/optimized_scatter.R`: R wrappers for the optimized C++ functions
- `R/PlotFlowFrame_optimized.R`: An optimized version of the PlotFlowFrame function

## Usage

To use the optimized functions, simply replace the original function calls with their optimized counterparts:

```r
# Original
rgbwt <- scatter_points_rgbwt(xy, out_size = c(512, 512), xlim = c(0, 1), ylim = c(0, 1))
rgba <- rgbwt_to_rgba_int(rgbwt)
raster <- rgba_int_to_raster(rgba)

# Optimized
raw_output <- scatter_points_optimized(xy, out_size = c(512, 512), xlim = c(0, 1), ylim = c(0, 1))
raster <- raw_to_raster(raw_output, c(512, 512))
```

For the main plotting function, use `PlotFlowFrame_optimized` instead of `PlotFlowFrame`:

```r
# Original
PlotFlowFrame(obj, channels, ...)

# Optimized
PlotFlowFrame_optimized(obj, channels, ..., use_raw_arrays = TRUE)
```

The `use_raw_arrays` parameter controls whether to use the optimized functions or fall back to the original implementation. This can be useful for debugging or comparison.

## Performance Comparison

The optimized functions can be significantly faster than the original implementation, especially for large datasets or high-resolution plots. See the `examples/optimized_example.R` script for a performance comparison.

## Limitations

The current optimization has some limitations:

1. Not all features of the original implementation are fully optimized. For complex operations like multi-population plotting with blending, the optimized version may fall back to the original implementation for some steps.
2. The optimized functions may have slightly different behavior in edge cases, although they should produce visually identical results in most cases.

## Future Improvements

Future improvements could include:

1. Implementing SIMD (Single Instruction Multiple Data) instructions for even better performance
2. Optimizing the memory layout for better cache locality
3. Extending the parallel processing to more operations
4. Adding GPU acceleration for very large datasets
