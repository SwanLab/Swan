#include "mex.h"
#include "math.h"
#include "func_cross.c"
#include "func_BarycentricCoordinatesTetrahedron.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    int odims[2]={1, 4};
    
    double *Lambda;
    double *p1, *p2, *p3, *p4, *RX, *RY, *RZ;
    
    p1=(double *)mxGetData(prhs[0]);
    p2=(double *)mxGetData(prhs[1]);
    p3=(double *)mxGetData(prhs[2]);
    p4=(double *)mxGetData(prhs[3]);
    RX=(double *)mxGetData(prhs[4]);
    RY=(double *)mxGetData(prhs[5]);
    RZ=(double *)mxGetData(prhs[6]);
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    Lambda=(double *)mxGetData(plhs[0]);
    
    BarycentricCoordinatesTetrahedron(p1, p2, p3, p4, RX[0], RY[0], RZ[0], Lambda);
}
