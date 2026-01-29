# Example script demonstrating the use of GPU-accelerated functions
# This script compares the performance of the original, optimized, and GPU functions

# Load required libraries
library(speedyflowplot)
library(flowCore)
library(microbenchmark)

# Create a simple example dataset
set.seed(42)
n_points <- 100000
xy <- matrix(runif(n_points * 2, 0, 1), ncol = 2)
colnames(xy) <- c("FSC-A", "SSC-A")

# Create a flowFrame
ff <- flowFrame(xy)

# Define colors for points
colors <- rainbow(10)
pop_indices <- list(
  "pop1" = sample(1:n_points, n_points * 0.3),
  "pop2" = sample(1:n_points, n_points * 0.2),
  "pop3" = sample(1:n_points, n_points * 0.1)
)

# Define population colors
pop_colors <- list(
  "pop1" = "#FF0000",
  "pop2" = "#00FF00",
  "pop3" = "#0000FF"
)

# Check if GPU is available
cat("Checking if GPU is available...\n")
gpu_available <- gpu_is_available()
cat("GPU available:", gpu_available, "\n\n")

# Compare performance of scatter_points functions
cat("Comparing scatter_points functions...\n")
benchmark_scatter <- microbenchmark(
  original = {
    rgbwt <- scatter_points_rgbwt(
      xy = xy,
      out_size = c(512, 512),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    rgba <- rgbwt_to_rgba_int(rgbwt)
    raster <- rgba_int_to_raster(rgba)
  },
  optimized = {
    raw_output <- scatter_points_optimized(
      xy = xy,
      out_size = c(512, 512),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    raster <- raw_to_raster(raw_output, c(512, 512))
  },
  gpu = {
    raw_output <- scatter_points_gpu(
      xy = xy,
      out_size = c(512, 512),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    raster <- raw_to_raster(raw_output, c(512, 512))
  },
  times = 5
)

print(benchmark_scatter)
cat("\n")

# Compare performance of PlotFlowFrame functions
cat("Comparing PlotFlowFrame functions...\n")
benchmark_plot <- microbenchmark(
  original = {
    pdf(NULL) # Redirect plot output to null device
    PlotFlowFrame(
      obj = ff,
      channels = c("FSC-A", "SSC-A"),
      pop_idex.list = pop_indices,
      population_colors = pop_colors,
      res = 256
    )
    dev.off()
  },
  optimized = {
    pdf(NULL) # Redirect plot output to null device
    PlotFlowFrame_optimized(
      obj = ff,
      channels = c("FSC-A", "SSC-A"),
      pop_idex.list = pop_indices,
      population_colors = pop_colors,
      res = 256,
      use_raw_arrays = TRUE
    )
    dev.off()
  },
  gpu = {
    pdf(NULL) # Redirect plot output to null device
    PlotFlowFrame_gpu(
      obj = ff,
      channels = c("FSC-A", "SSC-A"),
      pop_idex.list = pop_indices,
      population_colors = pop_colors,
      res = 256,
      use_raw_arrays = TRUE,
      use_gpu = TRUE
    )
    dev.off()
  },
  times = 3
)

print(benchmark_plot)
cat("\n")

# Create a visual example
par(mfrow = c(1, 3))

# Original version
PlotFlowFrame(
  obj = ff,
  channels = c("FSC-A", "SSC-A"),
  pop_idex.list = pop_indices,
  population_colors = pop_colors,
  title = "Original",
  res = 256
)

# Optimized version
PlotFlowFrame_optimized(
  obj = ff,
  channels = c("FSC-A", "SSC-A"),
  pop_idex.list = pop_indices,
  population_colors = pop_colors,
  title = "Optimized",
  res = 256,
  use_raw_arrays = TRUE
)

# GPU version
PlotFlowFrame_gpu(
  obj = ff,
  channels = c("FSC-A", "SSC-A"),
  pop_idex.list = pop_indices,
  population_colors = pop_colors,
  title = "GPU (Placeholder)",
  res = 256,
  use_raw_arrays = TRUE,
  use_gpu = TRUE
)

# Reset plot layout
par(mfrow = c(1, 1))

cat("Note: The GPU version is currently a placeholder that falls back to the optimized CPU implementation.\n")
cat("In a full implementation, the GPU version would use OpenCL or CUDA for GPU acceleration.\n")
cat("See docs/gpu_acceleration.md for details on how GPU acceleration would be implemented.\n")
