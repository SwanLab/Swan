#include "mex.h"
#include "math.h"
#include "func_BarycentricCoordinatesTriangle.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{
    /*  Size of output */
    int odims[2]={1, 3};
    
    double *Lambda;
    double *VN0, *VN1, *VN2,*RX, *RY;
            
    VN0=(double *)mxGetData(prhs[0]);
    VN1=(double *)mxGetData(prhs[1]);
    VN2=(double *)mxGetData(prhs[2]);
    RX=(double *)mxGetData(prhs[3]);
    RY=(double *)mxGetData(prhs[4]);
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    Lambda=(double *)mxGetData(plhs[0]);
    
    BarycentricCoordinatesTriangle(VN0,VN1,VN2,RX[0],RY[0],Lambda);
}
