#include "mex.h"
// #include <omp.h>
#include <stdio.h>



void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex  -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims DiagOperate.c
    mwIndex k,global_idx;
// // 
    mwIndex *C_S = mxGetIr(prhs[0]);
    mwIndex *starts_S = mxGetJc(prhs[0]);
    double  *V_S = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    double* vec = mxGetPr(prhs[1]);
    
    char* op = *((char*)mxGetData(prhs[2]));
    mwIndex m = mxGetN(prhs[1]);
    if (m==1){
        m = mxGetM(prhs[1]);
    }
    
    
    if (op=='A'){ //ADD
        if (m!=n){
            printf("DiagOperate: Matrix and vector not in the same size - doing nothing!!!%d,%d",n,m);
            return;
        }
        for (k = 0 ; k < n ; k++){
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                if (C_S[global_idx] == k){
                    V_S[global_idx] += vec[k];
                    break;
                }
            }
        }
    }else if (op=='E'){ //Exchage
        if (m!=n){
            printf("DiagOperate: Matrix and vector not in the same size - doing nothing!!!%d,%d",n,m);
            return;
        }
        for (k = 0 ; k < n ; k++){   
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                if (C_S[global_idx] == k){
                    V_S[global_idx] = vec[k];
                    break;
                }
            }
        }
    }else if (op=='R'){ // multiply from right
        if (m!=mxGetM(prhs[0])){
            printf("DiagOperate(R): Matrix and vector not in the same size - doing nothing!!!%d,%d",mxGetM(prhs[0]),m);
            return;
        }
        for (k = 0 ; k < n ; k++){   
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                V_S[global_idx] *= vec[C_S[global_idx]];
            }
        }
    
    }else if (op=='L'){ // multiply from left
        if (m!=n){
            printf("DiagOperate(L): Matrix and vector not in the same size - doing nothing!!!%d,%d",n,m);
            return;
        }
        for (k = 0 ; k < n ; k++){   
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
                V_S[global_idx] *= vec[k];
            }
        }
    
    }
    
}
 