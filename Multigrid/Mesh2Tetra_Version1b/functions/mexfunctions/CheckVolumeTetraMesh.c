#include "mex.h"
#include "math.h"
#include "func_cross.c"
#define absd(a)	          (((a) < 0)   ? -(a) : (a))

double CheckVolumeTetraMesh(double *V, int *sizeV, double *F, int *sizeF)
{
    int i;
    double volume=0;
    double k[3],a[3],b[3],c[3],d[3];
    int Fi1, Fi2, Fi3,Fi4;
    double v;
    for (i=0; i<sizeF[0]; i++)
    {
        Fi1=(int)F[i]-1;
        Fi2=(int)F[i+sizeF[0]]-1;
        Fi3=(int)F[i+2*sizeF[0]]-1;
        Fi4=(int)F[i+3*sizeF[0]]-1;
        
        a[0]=V[Fi1];  a[1]=V[Fi1+sizeV[0]];  a[2]=V[Fi1+2*sizeV[0]]; 
        b[0]=V[Fi2];  b[1]=V[Fi2+sizeV[0]];  b[2]=V[Fi2+2*sizeV[0]]; 
        c[0]=V[Fi3];  c[1]=V[Fi3+sizeV[0]];  c[2]=V[Fi3+2*sizeV[0]]; 
        d[0]=V[Fi4];  d[1]=V[Fi4+sizeV[0]];  d[2]=V[Fi4+2*sizeV[0]]; 

        a[0]-=d[0]; a[1]-=d[1]; a[2]-=d[2];
        b[0]-=d[0]; b[1]-=d[1]; b[2]-=d[2];
        c[0]-=d[0]; c[1]-=d[1]; c[2]-=d[2];
        
        cross(b, c, k);
        v = absd(a[0]*k[0]+a[1]*k[1]+a[2]*k[2])/6;
        volume+=v;
    }
    return volume;
}

        
     
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*  Size of output */
    double *V,*F, *Volume;
    int odims[2]={1, 1};
    int sizeV[3];
    const mwSize *sizeVc;
    int sizeF[3];
    const mwSize *sizeFc;
        

    V=(double *)mxGetData(prhs[0]);
    F=(double *)mxGetData(prhs[1]);
    sizeVc= mxGetDimensions(prhs[0]); 
    sizeFc= mxGetDimensions(prhs[1]);
    sizeV[0]=sizeVc[0]; sizeV[1]=sizeVc[1]; sizeV[2]=sizeVc[2];
    sizeF[0]=sizeFc[0]; sizeF[1]=sizeFc[1]; sizeF[2]=sizeFc[2];
    
    /* Create output array */
    plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
    Volume=(double *)mxGetData(plhs[0]);
    Volume[0]=CheckVolumeTetraMesh(V, sizeV, F, sizeF);
}
 

