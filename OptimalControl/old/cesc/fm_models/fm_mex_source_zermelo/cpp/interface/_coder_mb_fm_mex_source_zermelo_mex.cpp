//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_mb_fm_mex_source_zermelo_mex.cpp
//
// Code generation for function 'mb_fm_mex_source_zermelo'
//

// Include files
#include "_coder_mb_fm_mex_source_zermelo_mex.h"
#include "_coder_mb_fm_mex_source_zermelo_api.h"

// Function Definitions
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&mb_fm_mex_source_zermelo_atexit);
  // Module initialization.
  mb_fm_mex_source_zermelo_initialize();
  // Dispatch the entry-point.
  unsafe_mb_fm_mex_source_zermelo_mexFunction(nlhs, plhs, nrhs, prhs);
  // Module termination.
  mb_fm_mex_source_zermelo_terminate();
}

emlrtCTX mexFunctionCreateRootTLS()
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, nullptr, 1,
                           nullptr, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

void unsafe_mb_fm_mex_source_zermelo_mexFunction(int32_T nlhs, mxArray *plhs[2],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[2])
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  const mxArray *b_prhs[2];
  const mxArray *outputs[2];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  // Check for proper number of arguments.
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        24, "mb_fm_mex_source_zermelo");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "mb_fm_mex_source_zermelo");
  }
  // Call the function.
  b_prhs[0] = prhs[0];
  b_prhs[1] = prhs[1];
  mb_fm_mex_source_zermelo_api(b_prhs, nlhs, outputs);
  // Copy over outputs to the caller.
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

// End of code generation (_coder_mb_fm_mex_source_zermelo_mex.cpp)
