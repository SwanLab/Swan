#include "mex.h"
#include "math.h"
#define mind(a, b)        (((a) < (b)) ?  (a) : (b))
#define maxd(a, b)        (((a) > (b)) ?  (a) : (b))
#define absd(a)	          (((a) < 0)   ? -(a) : (a))
#include "func_BarycentricCoordinatesTriangle.c"
#include "func_CheckInsideFace.c"
#include "func_LineTriangleIntersection.c"
#include "func_TriangleTriangleIntersection.c"
            
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    int odims[2]={1, 1};
    
    double *P1, *P2, *P3, *O1, *O2, *O3;
    bool *IgnoreCorners;
    double *interO;
    
    P1=(double *)mxGetData(prhs[0]);
    P2=(double *)mxGetData(prhs[1]);
    P3=(double *)mxGetData(prhs[2]);
    O1=(double *)mxGetData(prhs[3]);
    O2=(double *)mxGetData(prhs[4]);
    O3=(double *)mxGetData(prhs[5]);
    IgnoreCorners=(bool *)mxGetData(prhs[6]);
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    interO=(double *)mxGetData(plhs[0]);
    if(TriangleTriangleIntersection(P1,P2,P3,O1,O2,O3,IgnoreCorners[0])) {
        interO[0]=1.0;
    }
    else {
        interO[0]=0.0;
    }
}
