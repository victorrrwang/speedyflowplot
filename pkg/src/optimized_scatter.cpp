/*
 * This file is part of speedyflowplot.
 *
 * Copyright (C) 2024 Victor Wang <VictorrWang@gmail.com>
 *
 * speedyflowplot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * speedyflowplot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with speedyflowplot. If not, see <https://www.gnu.org/licenses/>.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <cmath>

#include "macros.h"
#include "scatters.h"

// Optimized version that works directly with raw uint8_t arrays
// This avoids unnecessary data copies between R and C
extern "C" SEXP scatter_points_raw(SEXP xy_r, SEXP xlim_r, SEXP ylim_r, 
                                  SEXP out_size_r, SEXP rgba_r, SEXP raw_output_r) {
    // Extract dimensions
    int size_out_x = INTEGER(out_size_r)[0];
    int size_out_y = INTEGER(out_size_r)[1];
    int size_out = size_out_x * size_out_y;
    
    // Get xy data
    double *xy = REAL(xy_r);
    int n_points = LENGTH(xy_r) / 2;  // xy has 2 columns
    
    // Get limits
    double *xlim = REAL(xlim_r);
    double *ylim = REAL(ylim_r);
    
    // Get raw output buffer
    uint8_t *output = RAW(raw_output_r);
    
    // Check if output buffer has the right size (RGBA for each pixel)
    if (LENGTH(raw_output_r) != size_out * 4) {
        error("Output buffer size mismatch");
    }
    
    // Initialize output to transparent black
    std::fill(output, output + size_out * 4, 0);
    
    // Calculate bin sizes for x and y
    double x_begin = xlim[0];
    double x_end = xlim[1];
    double x_bin = (size_out_x - 1) / (x_end - x_begin);
    
    double y_begin = ylim[1];  // Note: y is flipped in the original code
    double y_end = ylim[0];
    double y_bin = (size_out_y - 1) / (y_end - y_begin);
    
    // Get RGBA color
    double *rgba = REAL(rgba_r);
    bool single_color = (LENGTH(rgba_r) == 4);
    
    // Process points
    for (int i = 0; i < n_points; ++i) {
        // Calculate pixel coordinates
        int x = (int)((xy[i] - x_begin) * x_bin);
        int y = (int)((xy[i + n_points] - y_begin) * y_bin);
        
        // Skip if outside the bitmap
        if (x < 0 || x >= size_out_x || y < 0 || y >= size_out_y)
            continue;
        
        // Calculate pixel offset
        int offset = (y * size_out_x + x) * 4;
        
        // Get color for this point
        double r, g, b, a;
        if (single_color) {
            r = rgba[0];
            g = rgba[1];
            b = rgba[2];
            a = rgba[3];
        } else {
            r = rgba[i * 4 + 0];
            g = rgba[i * 4 + 1];
            b = rgba[i * 4 + 2];
            a = rgba[i * 4 + 3];
        }
        
        // Alpha blending with existing pixel
        double existing_a = output[offset + 3] / 255.0;
        double new_a = a / 255.0;
        double result_a = existing_a + new_a * (1.0 - existing_a);
        
        if (result_a > 0) {
            // Blend colors
            output[offset + 0] = (uint8_t)((output[offset + 0] * existing_a + r * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 1] = (uint8_t)((output[offset + 1] * existing_a + g * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 2] = (uint8_t)((output[offset + 2] * existing_a + b * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 3] = (uint8_t)(result_a * 255.0);
        }
    }
    
    return R_NilValue;
}

// Optimized version for indexed colors
extern "C" SEXP scatter_indexed_raw(SEXP xy_r, SEXP xlim_r, SEXP ylim_r, 
                                   SEXP out_size_r, SEXP palette_r, SEXP map_r, SEXP raw_output_r) {
    // Extract dimensions
    int size_out_x = INTEGER(out_size_r)[0];
    int size_out_y = INTEGER(out_size_r)[1];
    int size_out = size_out_x * size_out_y;
    
    // Get xy data
    double *xy = REAL(xy_r);
    int n_points = LENGTH(xy_r) / 2;  // xy has 2 columns
    
    // Get limits
    double *xlim = REAL(xlim_r);
    double *ylim = REAL(ylim_r);
    
    // Get raw output buffer
    uint8_t *output = RAW(raw_output_r);
    
    // Check if output buffer has the right size (RGBA for each pixel)
    if (LENGTH(raw_output_r) != size_out * 4) {
        error("Output buffer size mismatch");
    }
    
    // Initialize output to transparent black
    std::fill(output, output + size_out * 4, 0);
    
    // Calculate bin sizes for x and y
    double x_begin = xlim[0];
    double x_end = xlim[1];
    double x_bin = (size_out_x - 1) / (x_end - x_begin);
    
    double y_begin = ylim[1];  // Note: y is flipped in the original code
    double y_end = ylim[0];
    double y_bin = (size_out_y - 1) / (y_end - y_begin);
    
    // Get palette and map
    double *palette = REAL(palette_r);
    int *map = INTEGER(map_r);
    int palette_size = LENGTH(palette_r) / 4;  // Each color has 4 components
    
    // Process points
    for (int i = 0; i < n_points; ++i) {
        // Calculate pixel coordinates
        int x = (int)((xy[i] - x_begin) * x_bin);
        int y = (int)((xy[i + n_points] - y_begin) * y_bin);
        
        // Skip if outside the bitmap
        if (x < 0 || x >= size_out_x || y < 0 || y >= size_out_y)
            continue;
        
        // Calculate pixel offset
        int offset = (y * size_out_x + x) * 4;
        
        // Get color index for this point
        int color_idx = map[i] - 1;  // R indices are 1-based
        
        // Skip if index is out of range
        if (color_idx < 0 || color_idx >= palette_size)
            continue;
        
        // Get color from palette
        double r = palette[color_idx * 4 + 0];
        double g = palette[color_idx * 4 + 1];
        double b = palette[color_idx * 4 + 2];
        double a = palette[color_idx * 4 + 3];
        
        // Alpha blending with existing pixel
        double existing_a = output[offset + 3] / 255.0;
        double new_a = a / 255.0;
        double result_a = existing_a + new_a * (1.0 - existing_a);
        
        if (result_a > 0) {
            // Blend colors
            output[offset + 0] = (uint8_t)((output[offset + 0] * existing_a + r * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 1] = (uint8_t)((output[offset + 1] * existing_a + g * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 2] = (uint8_t)((output[offset + 2] * existing_a + b * new_a * (1.0 - existing_a)) / result_a);
            output[offset + 3] = (uint8_t)(result_a * 255.0);
        }
    }
    
    return R_NilValue;
}

// Apply kernel to raw RGBA bitmap
extern "C" SEXP apply_kernel_raw(SEXP raw_input_r, SEXP out_size_r, SEXP kernel_r, SEXP raw_output_r) {
    // Extract dimensions
    int size_out_x = INTEGER(out_size_r)[0];
    int size_out_y = INTEGER(out_size_r)[1];
    int size_out = size_out_x * size_out_y;
    
    // Get input and output buffers
    uint8_t *input = RAW(raw_input_r);
    uint8_t *output = RAW(raw_output_r);
    
    // Check if buffers have the right size
    if (LENGTH(raw_input_r) != size_out * 4 || LENGTH(raw_output_r) != size_out * 4) {
        error("Buffer size mismatch");
    }
    
    // Get kernel
    double *kernel = REAL(kernel_r);
    int kernel_size = (int)sqrt(LENGTH(kernel_r));
    int radius = kernel_size / 2;
    
    // Initialize output to transparent black
    std::fill(output, output + size_out * 4, 0);
    
    // Apply kernel to each pixel
    for (int y = 0; y < size_out_y; ++y) {
        for (int x = 0; x < size_out_x; ++x) {
            double r = 0, g = 0, b = 0, a = 0;
            
            // Apply kernel
            for (int ky = -radius; ky <= radius; ++ky) {
                for (int kx = -radius; kx <= radius; ++kx) {
                    int px = x + kx;
                    int py = y + ky;
                    
                    // Skip if outside the bitmap
                    if (px < 0 || px >= size_out_x || py < 0 || py >= size_out_y)
                        continue;
                    
                    // Get kernel weight
                    double weight = kernel[(ky + radius) * kernel_size + (kx + radius)];
                    
                    // Get pixel color
                    int offset = (py * size_out_x + px) * 4;
                    r += input[offset + 0] * weight;
                    g += input[offset + 1] * weight;
                    b += input[offset + 2] * weight;
                    a += input[offset + 3] * weight;
                }
            }
            
            // Write result
            int offset = (y * size_out_x + x) * 4;
            output[offset + 0] = (uint8_t)std::min(255.0, std::max(0.0, r));
            output[offset + 1] = (uint8_t)std::min(255.0, std::max(0.0, g));
            output[offset + 2] = (uint8_t)std::min(255.0, std::max(0.0, b));
            output[offset + 3] = (uint8_t)std::min(255.0, std::max(0.0, a));
        }
    }
    
    return R_NilValue;
}

// Register the new functions
static const R_CallMethodDef callMethods[] = {
    {"scatter_points_raw", (DL_FUNC) &scatter_points_raw, 6},
    {"scatter_indexed_raw", (DL_FUNC) &scatter_indexed_raw, 7},
    {"apply_kernel_raw", (DL_FUNC) &apply_kernel_raw, 4},
    {NULL, NULL, 0}
};

extern "C" void R_init_optimized_scatter(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
