//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_mb_fm_mex_source_zermelo_info.cpp
//
// Code generation for function 'mb_fm_mex_source_zermelo'
//

// Include files
#include "_coder_mb_fm_mex_source_zermelo_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

// Function Declarations
static const mxArray *emlrtMexFcnResolvedFunctionsInfo();

// Function Definitions
static const mxArray *emlrtMexFcnResolvedFunctionsInfo()
{
  const mxArray *nameCaptureInfo;
  const char_T *data[5]{
      "789cc554cb4ec240141d0c3e36282bbfc19569e55170a5a13e13450326268e29657aabd5"
      "ce0ce9b408ee4ddcb9d2cf72e1d6855b3fc1a5522850e2041222decd"
      "9d93939973ee99c9a0c4c1510221b48cba955eecf654847b7d0ec56b944f487a54f32819"
      "db17f14fbd4e38f3a1e577013329f4775a9c3acc647eb5dd00e481e0",
      "6e13ac90b11d17aa0e85ca3038ee20ba3b44f54187eaac4bd7406e2b0145deb518387487"
      "413f8f4fc9bcc909f37890e4911ee12f762e4b9bf84c802770cb6c3a"
      "0d5c66a07b4e13b0ce494081f95ce03dc7df0feab87267325c6ef80e35ddd24f701e7731"
      "0141b04d0dca2d7045b882962178e01130eec1a3e0726c39b68d69dd",
      "f8955da743733f4e39f7da98b9235eb485912916c056f2d9a2aa65b24a319f2b281a5188"
      "aadaaaa98169193726897cd5a6f4b520f5d5652c1ed45d18e4f035a5"
      "deb3542fcecff0fe278d3cf61e6a92395726cc41f62fa4d052d85f5edf436a567a871f5b"
      "6fb3d48beabff45a92f3267dc7ab12bdf408af072ca7f1ea39dbe679",
      "ce3772ea15390df4818f93313ae37c2009feebf3bf012de1b6d5", ""};
  nameCaptureInfo = nullptr;
  emlrtNameCaptureMxArrayR2016a(&data[0], 1832U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties()
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[7]{
      "Version",      "ResolvedFunctions", "Checksum",    "EntryPoints",
      "CoverageInfo", "IsPolymorphic",     "PropertyList"};
  const char_T *epFieldName[6]{
      "Name",           "NumberOfInputs", "NumberOfOutputs",
      "ConstantInputs", "FullPath",       "TimeStamp"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 6, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("mb_fm_mex_source_zermelo"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString(
          "C:"
          "\\Users\\xavip\\OneDrive\\Documentos\\GitHub\\Swan\\OptimalControl\\"
          "cesc\\fm_models\\fm_mex_source_zermelo\\diff\\mb_fm_mex_so"
          "urce_zermelo.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739764.06609953707));
  xResult =
      emlrtCreateStructMatrix(1, 1, 7, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("9.14.0.2891782 (R2023a) Update 8"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("P2BaLQnArfOPfnWYWVXwt"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

// End of code generation (_coder_mb_fm_mex_source_zermelo_info.cpp)
