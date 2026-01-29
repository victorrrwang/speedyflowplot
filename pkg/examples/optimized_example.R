# Example script demonstrating the use of optimized functions
# This script compares the performance of the original and optimized functions

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

# Compare performance of scatter_points_rgbwt vs scatter_points_optimized
cat("Comparing scatter_points_rgbwt vs scatter_points_optimized...\n")
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
  times = 10
)

print(benchmark_scatter)
cat("\n")

# Compare performance of PlotFlowFrame vs PlotFlowFrame_optimized
cat("Comparing PlotFlowFrame vs PlotFlowFrame_optimized...\n")
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
  times = 5
)

print(benchmark_plot)
cat("\n")

# Create a visual example
par(mfrow = c(1, 2))

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

# Reset plot layout
par(mfrow = c(1, 1))

cat("The optimized version should be significantly faster while producing identical visual results.\n")
cat("The performance improvement comes from avoiding unnecessary data copies between R and C,\n")
cat("and working directly with raw uint8_t arrays instead of converting between different formats.\n")
