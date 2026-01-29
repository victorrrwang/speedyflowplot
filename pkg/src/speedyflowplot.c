#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Function declarations for .C interface
void scatter_indexed_rgbwt(int *dimen, float *xlim, float *ylim, float *palette, float *fRGBWT, int *map, float *xy);
void scatter_singlecolor_rgbwt(int *dimen, float *xlim, float *ylim, float *RGBA, float *fRGBWT, float *xy);
void scatter_multicolor_rgbwt(int *dimen, float *xlim, float *ylim, float *RGBA, float *fRGBWT, float *xy);
void kernel_rgbwt(int *dimen, float *kernel, float *blurred_fRGBWT, float *fRGBWT);

// Function declarations for .Call interface
SEXP scatter_points_raw(SEXP xy_r, SEXP xlim_r, SEXP ylim_r, SEXP out_size_r, SEXP rgba_r, SEXP raw_output_r);
SEXP scatter_indexed_raw(SEXP xy_r, SEXP xlim_r, SEXP ylim_r, SEXP out_size_r, SEXP palette_r, SEXP map_r, SEXP raw_output_r);
SEXP apply_kernel_raw(SEXP raw_input_r, SEXP out_size_r, SEXP kernel_r, SEXP raw_output_r);

// Register .C interface routines
static const R_CMethodDef cMethods[] = {
    {"scatter_indexed_rgbwt", (DL_FUNC) &scatter_indexed_rgbwt, 7},
    {"scatter_singlecolor_rgbwt", (DL_FUNC) &scatter_singlecolor_rgbwt, 6},
    {"scatter_multicolor_rgbwt", (DL_FUNC) &scatter_multicolor_rgbwt, 6},
    {"kernel_rgbwt", (DL_FUNC) &kernel_rgbwt, 4},
    {NULL, NULL, 0}
};

// Register .Call interface routines
static const R_CallMethodDef callMethods[] = {
    {"scatter_points_raw", (DL_FUNC) &scatter_points_raw, 6},
    {"scatter_indexed_raw", (DL_FUNC) &scatter_indexed_raw, 7},
    {"apply_kernel_raw", (DL_FUNC) &apply_kernel_raw, 4},
    {NULL, NULL, 0}
};

void R_init_speedyflowplot(DllInfo *info) {
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
