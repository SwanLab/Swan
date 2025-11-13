//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// mb_fm_mex_source_zermelo.h
//
// Code generation for function 'mb_fm_mex_source_zermelo'
//

#ifndef MB_FM_MEX_SOURCE_ZERMELO_H
#define MB_FM_MEX_SOURCE_ZERMELO_H

// Include files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Custom Header Code
#include <cstdint>
#define CHAR16_T uint16_t
#include "mex.h"
// Function Declarations
extern void mb_fm_mex_source_zermelo(const double states[2], double controls,
                                     double statesdot[2],
                                     double j_statesdot[6]);

extern void mb_fm_mex_source_zermelo_initialize();

extern void mb_fm_mex_source_zermelo_terminate();

#endif
// End of code generation (mb_fm_mex_source_zermelo.h)
