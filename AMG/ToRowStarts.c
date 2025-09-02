/*---------------------------------------------------------------------------
Copyright (2010): Eran Treister and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/
#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int i;
    double* In = mxGetPr(prhs[0]);
    int* Out = 0;
    int n = (int)*mxGetPr(prhs[1]);
    int nnz = (int)*mxGetPr(prhs[2]);
    const int dims[] = {1,n}; 
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    /*Get a pointer to the data space in our newly allocated memory*/
    Out = (int*)mxGetData(plhs[0]);
    for (i = 0 ; i < n ; i++){
        Out[i] = -1;
    }
    Out[(int)In[0]-1] = 0;
    for ( i = 1 ; i < nnz ; i++){
        if (In[i] > In[i-1] + 1e-12){
            Out[(int)In[i]-1] = i;  
        }
	}    
}