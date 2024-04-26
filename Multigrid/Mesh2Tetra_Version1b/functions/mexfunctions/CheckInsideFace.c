#include "mex.h"
#include "math.h"
#include "func_BarycentricCoordinatesTriangle.c"
#include "func_CheckInsideFace.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    int odims[2]={1, 1};
    
    double *VN0, *VN1, *VN2, *RX, *RY;
    double *MI, *MA;
    double *interO;
    bool *IgnoreCorners;
    
    VN0=(double *)mxGetData(prhs[0]);
    VN1=(double *)mxGetData(prhs[1]);
    VN2=(double *)mxGetData(prhs[2]);
    RX=(double *)mxGetData(prhs[3]);
    RY=(double *)mxGetData(prhs[4]);
    MI=(double *)mxGetData(prhs[5]);
    MA=(double *)mxGetData(prhs[6]);
    IgnoreCorners=(bool *)mxGetData(prhs[7]);
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    interO=(double *)mxGetData(plhs[0]);
    if(CheckInsideFace(VN0, VN1, VN2, RX[0], RY[0], MI, MA,  IgnoreCorners[0])) {
        interO[0]=1.0;
    }
    else {
        interO[0]=0.0;
    }
 }
