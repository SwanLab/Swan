#include "mex.h"
#include <omp.h>
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex -largeArrayDims MulticolorGaussSeidelSMatTranspose.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"
    mwIndex i,j,i_color;
    int id, Nthrds, istart, iend,k;
    long idx_in_color;
    mwIndex *C_t = mxGetIr(prhs[0]);
    mwIndex *starts_t = mxGetJc(prhs[0]);
    double* vals_t = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    double temp;
    
    double* x_in = mxGetPr(prhs[1]);
    double* b = mxGetPr(prhs[2]);
    double* invDiag = mxGetPr(prhs[3]);
    int nu = (int)(*mxGetPr(prhs[4]));
    double* indicesOrder = mxGetPr(prhs[5]);
    double* colorStarts = mxGetPr(prhs[6]);
    mwIndex n_colors = max(mxGetM(prhs[6]),mxGetN(prhs[6]))-1;
    /* Output Variables */
    double* x_out = 0;
    double* aux = (double*)malloc(n*sizeof(double));
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_out = mxGetPr(plhs[0]);
    /* Program */
    
    #pragma omp parallel shared(colorStarts, indicesOrder,vals_t,C_t,starts_t,invDiag,x_in,b,aux) private(i_color,idx_in_color,k,i,j,temp,id, Nthrds, istart, iend) num_threads(omp_get_num_procs()/2)
    {
        //printf("%d",);
        id = omp_get_thread_num();
        Nthrds = omp_get_num_threads();
        istart = id * n / Nthrds;
        iend = (id+1) * n / Nthrds;
        if (id == Nthrds-1)iend = n;
        for ( i = istart ; i < iend ; ++i){
            aux[i] =  x_in[i];
        }
        for (k = 0 ; k < nu ; k++){
            #pragma omp barrier
            for ( i_color = 0 ; i_color < n_colors ; i_color++){
                istart = id * ((colorStarts[i_color+1] - colorStarts[i_color]) / Nthrds) + colorStarts[i_color];
                iend = (id+1) * ((colorStarts[i_color+1] - colorStarts[i_color]) / Nthrds) + colorStarts[i_color];
                if (id == Nthrds-1)iend = colorStarts[i_color+1];
                for ( idx_in_color = istart ; idx_in_color < iend ; ++idx_in_color){
                    i = indicesOrder[idx_in_color];
                    temp = 0.0;
                    for (j = starts_t[i] ; j < starts_t[i+1] ; ++j)
                    {
                        temp += vals_t[j]*aux[C_t[j]];
                    }
                    x_out[i] = aux[i] + (b[i] - temp)*invDiag[i];
                }
                #pragma omp barrier
                for ( idx_in_color = istart ; idx_in_color < iend ; ++idx_in_color){
                    i = indicesOrder[idx_in_color];
                    aux[i] = x_out[i];
                }
                #pragma omp barrier
           }
        }
    }
    free(aux);
}