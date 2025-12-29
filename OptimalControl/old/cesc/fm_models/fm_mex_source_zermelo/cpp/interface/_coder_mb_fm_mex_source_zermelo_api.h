//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_mb_fm_mex_source_zermelo_api.h
//
// Code generation for function 'mb_fm_mex_source_zermelo'
//

#ifndef _CODER_MB_FM_MEX_SOURCE_ZERMELO_API_H
#define _CODER_MB_FM_MEX_SOURCE_ZERMELO_API_H

// Include files
#include "emlrt.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void mb_fm_mex_source_zermelo(real_T states[2], real_T controls,
                              real_T statesdot[2], real_T j_statesdot[6]);

void mb_fm_mex_source_zermelo_api(const mxArray *const prhs[2], int32_T nlhs,
                                  const mxArray *plhs[2]);

void mb_fm_mex_source_zermelo_atexit();

void mb_fm_mex_source_zermelo_initialize();

void mb_fm_mex_source_zermelo_terminate();

void mb_fm_mex_source_zermelo_xil_shutdown();

void mb_fm_mex_source_zermelo_xil_terminate();

#endif
// End of code generation (_coder_mb_fm_mex_source_zermelo_api.h)
