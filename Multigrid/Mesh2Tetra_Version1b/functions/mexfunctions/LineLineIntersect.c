#include "mex.h"
#include "math.h"
#define mind(a, b)        (((a) < (b)) ?  (a) : (b))
#define maxd(a, b)        (((a) > (b)) ?  (a) : (b))
#define absd(a)	          (((a) < 0)   ? -(a) : (a))
#include "func_LineLineIntersect.c"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    int odims[2]={1, 1};
    
    double *P1, *P2, *O1, *O2;
    double *Inter, *Inter_x, *Inter_y;
    double Val[3];
    double P1a[2], P2a[2], O1a[2], O2a[2];
    
    P1=(double *)mxGetData(prhs[0]);
    P2=(double *)mxGetData(prhs[1]);
    O1=(double *)mxGetData(prhs[2]);
    O2=(double *)mxGetData(prhs[3]);
   
    
    P1a[0]=P1[0]; P1a[1]=P1[1];
    P2a[0]=P2[0]; P2a[1]=P2[1];
    O1a[0]=O1[0]; O1a[1]=O1[1];
    O2a[0]=O2[0]; O2a[1]=O2[1];

    LineLineIntersect( P1a,P2a,O1a,O2a,Val);
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    Inter  =(double *)mxGetData(plhs[0]);
    Inter[0]  =Val[0];
    if(nrhs>1)
    {
        plhs[1] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
        Inter_x=(double *)mxGetData(plhs[1]);
        Inter_y=(double *)mxGetData(plhs[2]);
        Inter_x[0]=Val[1];
        Inter_y[0]=Val[2];
    }
}
 
