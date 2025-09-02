#include "mex.h"
#include <omp.h>
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex -largeArrayDims GaussSeidelSMatTranspose.c
    int i;
    int id, Nthrds, istart, iend;
    mwIndex j,k;
    double temp;
    double omega = 0.8;
    mwIndex *C_t = mxGetIr(prhs[0]);
    mwIndex *starts_t = mxGetJc(prhs[0]);
    double* vals_t = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    double* x_in = mxGetPr(prhs[1]);
    double* b = mxGetPr(prhs[2]);
    double* invDiag = mxGetPr(prhs[3]);
    int nu = (int)(*mxGetPr(prhs[4]));
    /* Output Variables */
    double* x_out = 0;
    double* aux = 0;
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_out = mxGetPr(plhs[0]);
    if (nu>1){
        aux = (double*)malloc(n*sizeof(double));
    }else{
        aux = x_in;
    }
    #pragma omp parallel shared(vals_t,C_t,starts_t,invDiag,x_in,b,omega,aux) private(i,j,k,temp,id, Nthrds, istart, iend) num_threads(2)
    {
        id = omp_get_thread_num();
        Nthrds = omp_get_num_threads();
        istart = id * n / Nthrds;
        iend = (id+1) * n / Nthrds;
        if (id == Nthrds-1)iend = n;
        for (k = 0 ; k < nu ; k++){
            #pragma omp barrier
            if (nu>1){
                for ( i = istart ; i < iend ; ++i){
                    aux[i] = x_out[i];
                }
            }
            #pragma omp barrier
            for ( i = istart ; i < iend ; ++i){
                temp = 0.0;
                for (j = starts_t[i] ; j < starts_t[i+1] ; ++j)
                {
                    temp += vals_t[j]*aux[C_t[j]];
                }
                x_out[i] = aux[i] + omega*(b[i] - temp)*invDiag[i];
            }
        }
    }
    if (nu>1){
        free(aux);
    }
}