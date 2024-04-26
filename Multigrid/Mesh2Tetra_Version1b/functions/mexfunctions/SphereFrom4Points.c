#include "mex.h"
#include "math.h"
#include "func_Determinant.c"
#include "func_SphereFrom4Points.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    int odims[2]={1, 1};
    double *PA, *PB, *PC, *PD;
    double *Cx, *Cy, *Cz, *L;
    
    PA=(double *)mxGetData(prhs[0]);
    PB=(double *)mxGetData(prhs[1]);
    PC=(double *)mxGetData(prhs[2]);
    PD=(double *)mxGetData(prhs[3]);
    
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    Cx =(double *)mxGetData(plhs[0]);
    Cy =(double *)mxGetData(plhs[1]);
    Cz =(double *)mxGetData(plhs[2]);
    L = (double *)mxGetData(plhs[3]);
    
    SphereFrom4Points(PA, PB, PC, PD, Cx, Cy, Cz, L);
}

