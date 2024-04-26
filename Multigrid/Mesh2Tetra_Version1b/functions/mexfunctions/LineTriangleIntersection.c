#include "mex.h"
#include "math.h"
#define mind(a, b)        (((a) < (b)) ?  (a) : (b))
#define maxd(a, b)        (((a) > (b)) ?  (a) : (b))
#include "func_BarycentricCoordinatesTriangle.c"
#include "func_CheckInsideFace.c"
#include "func_LineTriangleIntersection.c"



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    int odims[2]={1, 1};
    
    double *A, *B, *C, *P1, *P2;
    bool *IgnoreCorners;
    bool Ignore;
	double *interO;
	double inter_xyz[3];
	double *inter_x;
	double *inter_y;
	double *inter_z;
	
    
    A=(double *)mxGetData(prhs[0]);
    B=(double *)mxGetData(prhs[1]);
    C=(double *)mxGetData(prhs[2]);
    P1=(double *)mxGetData(prhs[3]);
    P2=(double *)mxGetData(prhs[4]);
	if(nrhs>5)
	{
	    IgnoreCorners=(bool *)mxGetData(prhs[5]);
		Ignore=IgnoreCorners[0];
	}
    else
	{
		Ignore=true;
	}
	
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    interO=(double *)mxGetData(plhs[0]);
    if(LineTriangleIntersection(A, B, C, P1, P2, Ignore,inter_xyz)) {
        interO[0]=1.0;
    }
    else {
        interO[0]=0.0;
    }
	if(nlhs>1)
	{
		plhs[1] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
		plhs[2] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
		plhs[3] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
		inter_x=(double *)mxGetData(plhs[1]);
		inter_y=(double *)mxGetData(plhs[2]);
		inter_z=(double *)mxGetData(plhs[3]);
		inter_x[0]=inter_xyz[0];
		inter_y[0]=inter_xyz[1];
		inter_z[0]=inter_xyz[2];
	}
}
