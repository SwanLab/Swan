#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex -largeArrayDims GaussSeidelSMatTranspose.c
    int i,j,k;
    double temp;
    mwIndex *C_t = mxGetIr(prhs[0]);
    mwIndex *starts_t = mxGetJc(prhs[0]);
    double* vals_t = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    double* x_in = mxGetPr(prhs[1]);
    double* b = mxGetPr(prhs[2]);
    double* invDiag = mxGetPr(prhs[3]);
    char direction = *((char*)mxGetData(prhs[4]));
    
    /* Output Variables */
    double* x_out = 0;
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_out = mxGetPr(plhs[0]);
    
    for (i = 0 ; i < n ; ++i){
        x_out[i] = x_in[i];
    }
   
    if ((direction=='F')||(direction=='S')){
        for ( i = 0 ; i < n ; ++i){
            temp = 0.0;
            for (j = starts_t[i] ; j < starts_t[i+1] ; ++j)
            {
                temp += vals_t[j]*x_out[C_t[j]];
            }
            x_out[i] += (b[i] - temp)*invDiag[i];
        }
    }
    if ((direction=='B')||(direction=='S')){
        for ( i = n-1 ; i >= 0 ; --i){     
            temp = 0.0;
            for (j = starts_t[i] ; j < starts_t[i+1] ; ++j)
            {
                temp += vals_t[j]*x_out[C_t[j]];
            }
            x_out[i] += (b[i] - temp)*invDiag[i];
        }      
    }
}