#include "mex.h"
#include <omp.h>
#include <math.h>


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex  -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims NeighborhoodAgg.c
    mwIndex k,global_idx;
// // 
    mwIndex *C_S = mxGetIr(prhs[0]);
    mwIndex *starts_S = mxGetJc(prhs[0]);
    double  *V_S = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    double* rowQ;
    double* colQ;
    double* used;
    double* numAggregated_ans = 0;
    double* numAggs_ans = 0;
    mwIndex numAggregated = 0;
    mwIndex numAggs = 0;
    mwIndex neighbor;
    int Neighbors_Aggregated_flag;
    
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(n, 1, mxREAL);
    numAggs_ans = mxGetPr(plhs[0]);
    numAggregated_ans = mxGetPr(plhs[1]);
    used = mxGetPr(plhs[2]);
    rowQ = mxGetPr(plhs[3]);
    colQ = mxGetPr(plhs[4]);
    Neighbors_Aggregated_flag = 0;
    for (k = 0 ; k < n ; k++){
        used[k] = 0;
    }
    for (k = 0 ; k < n ; k++){
        Neighbors_Aggregated_flag = 0;
        for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
            if (used[C_S[global_idx]] == 1){
                Neighbors_Aggregated_flag = 1;
                break;
            }
        }
        if (Neighbors_Aggregated_flag==0){
            numAggs++; // starts from 1...
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                neighbor = C_S[global_idx];
                used[neighbor] = 1;
                rowQ[numAggregated] = neighbor+1; // conversion to MATLAB indices
                colQ[numAggregated] = numAggs;
                numAggregated++;
            }
        }
    }
    *numAggregated_ans = numAggregated;
    *numAggs_ans = numAggs;
    
}
 