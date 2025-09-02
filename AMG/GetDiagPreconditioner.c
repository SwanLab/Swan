#include "mex.h"
// #include <omp.h>
#include <stdio.h>



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex  -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims GetDiagPreconditioner.c
    mwIndex k,global_idx;
// // 
    mwIndex *C_S = mxGetIr(prhs[0]);
    mwIndex *starts_S = mxGetJc(prhs[0]);
    double  *V_S = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    double diag = 0;
    char op = *((char*)mxGetData(prhs[1]));
    double* x_out = 0;
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_out = mxGetPr(plhs[0]);
    
    
    
    if (op=='1'){ // L1
        for (k = 0 ; k < n ; k++){
            x_out[k] = 0;
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                x_out[k]+=fabs(V_S[global_idx]);
            }
            if (x_out[k] > 0.0){
                x_out[k] = 1/x_out[k];
            }
        }
    }else if (op=='S'){ //SPAI
        for (k = 0 ; k < n ; k++){
            x_out[k] = 0;
            diag = 0;
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                x_out[k] += (V_S[global_idx]*V_S[global_idx]);
                diag += (C_S[global_idx] == k)*V_S[global_idx];
            }
            if (x_out[k] > 0.0){
                x_out[k] = diag/x_out[k];
            }
        }
    }
}
 