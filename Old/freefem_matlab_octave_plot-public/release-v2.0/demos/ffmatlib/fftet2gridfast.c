/* Interpolates from a 3D tetrahedral mesh to a 2D grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2018-11-13
 *
 *   This file is part of the ffmatlib which is hosted at
 *   https://github.com/samplemaker/freefem_matlab_octave_plot
 *
 *   [w1, ...] = fftet2gridfast (x, y, z, tx, ty, tz, tu1, ...)
 *
 *   fftet2grid computes the function values w1, w2, ... over a mesh grid defined
 *   by the arguments x, y from a set of functions u1, u2, ... with values
 *   given on a tetrahedral mesh tx, ty, tz. The values are computed using first order
 *   or second order approximating basis functions (P1 or P2 - Lagrangian Finite
 *   Elements). The function values w1, w2, ... are real if tu1, tu2, ... are real
 *   or complex if tu1, tu2, ... are complex. The mesh grid x, y, z can be cartesian
 *   or curved. fftet2grid returns NaNs when an interpolation point is outside the
 *   tetrahedral mesh. fftetgridfast.c is a runtime optimized mex implementation of
 *   the function fftet2grid.m.
 *
 *   This code is compatible with OCTAVE and MATLAB before R2018
 *   (old Complex API)
 *
 *   Octave users on Linux with gcc compile the mex file with the command:
 *
 *       mkoctfile --mex -Wall fftet2gridfast.c
 *
 *   Matlab users on Windows with microsoft visual studio compile the mex
 *   file with the command:
 *
 *       mex fftet2gridfast.c -largeArrayDims
 *
 * Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 */


#include <math.h>
#include "mex.h"

#define min1(a,b) ( ((a)<(b)) ? (a):(b) )
#define max1(a,b) ( ((a)>(b)) ? (a):(b) )

#define P1(Va,Vb,Vc,Vd,u1,u2,u3,u4) (u1*Va+u2*Vb+u3*Vc+u4*Vd)

#define P2(Va,Vb,Vc,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10) u1*Va*(2.0*Va-1.0)+                      \
                                                    u2*Vb*(2.0*Vb-1.0)+                      \
                                                    u3*Vc*(2.0*Vc-1.0)+                      \
                                                    u4*(1.0-Va-Vb-Vc)*(1.0-2.0*(Va+Vb+Vc))+  \
                                                    u5*4.0*Va*Vb+                            \
                                                    u6*4.0*Va*Vc+                            \
                                                    u7*4.0*Va*(1.0-Va-Vb-Vc)+                \
                                                    u8*4.0*Vb*Vc+                            \
                                                    u9*4.0*Vb*(1.0-Va-Vb-Vc)+                \
                                                    u10*4.0*Vc*(1.0-Va-Vb-Vc)

static inline
double max4(double i, double j, double k, double l){
  return max1(max1(i,j), max1(k,l));
}

static inline
double min4(double i, double j, double k, double l){
  return min1(min1(i,j), min1(k,l));
}

typedef struct{
  double i,j,k;
}Vector;

static inline
double dotProduct(Vector a, Vector b){
  return a.i*b.i+a.j*b.j+a.k*b.k;
}

static inline
Vector crossProduct(Vector a, Vector b){
  Vector c={a.j*b.k-a.k*b.j, a.k*b.i-a.i*b.k, a.i*b.j-a.j*b.i};
  return c;
}

static inline
Vector minus(Vector a, Vector b){
  Vector c={a.i-b.i, a.j-b.j, a.k-b.k};
  return c;
}

typedef enum {
  P0,
  P1,
  P2
} feElement;

void
fftet2gridfast (double *x, double *y, double *z, double *tx, double *ty,
                double *tz, double **tuRe, double **tuIm, double **wRe,
                double **wIm, mwSize nOuts, mwSize nTet, mwSize nx,
                mwSize ny, feElement elementType){

  double *invV0=(double *)mxMalloc(nTet*sizeof(double));
  double init=mxGetNaN( );

  /*Triangle areas */
  mwSize j=0;
  mwSize nNodes=4*nTet;
  for (mwSize i=0; i<nTet; i++) {
    Vector ap={tx[j], ty[j], tz[j]},
           bp={tx[j+1], ty[j+1], tz[j+1]},
           cp={tx[j+2], ty[j+2], tz[j+2]},
           dp={tx[j+3], ty[j+3], tz[j+3]};
    invV0[i]=fabs(1.0/(dotProduct(crossProduct(minus(cp,ap),
                                               minus(bp,ap)),
                                               minus(dp,ap))));
    j=j+4;
  }
  /*For all grid points of the meshgrid */
  for (mwSize mx=0; mx<nx; mx++){
    for (mwSize my=0; my<ny; my++){
      mwSize ofs=(mx*ny)+my;
      for (mwSize nArg=0; nArg<nOuts; nArg++){
        if (tuIm[nArg] != NULL){
          *(wIm[nArg]+ofs)=init;
          *(wRe[nArg]+ofs)=init;
        }
        else{
          *(wRe[nArg]+ofs)=init;
        }
      }
      mwSize i=0, j=0;
      bool doSearchTet=true;
      while (doSearchTet && (i<nNodes)){
        /*If the point (X,Y) is outside a square defined by the tetrahedron
          vertices, the point can not be within the tetrahedron */
        bool presel=((x[ofs] <= max4(tx[i],tx[i+1],tx[i+2],tx[i+3])) &&
                     (x[ofs] >= min4(tx[i],tx[i+1],tx[i+2],tx[i+3])) &&
                     (y[ofs] <= max4(ty[i],ty[i+1],ty[i+2],ty[i+3])) &&
                     (y[ofs] >= min4(ty[i],ty[i+1],ty[i+2],ty[i+3])));
        /*Potential candiate - calculate Barycentric Coordinates */
        if (presel) {
          Vector ap={tx[i], ty[i], tz[i]},
                 bp={tx[i+1], ty[i+1], tz[i+1]},
                 cp={tx[i+2], ty[i+2], tz[i+2]},
                 dp={tx[i+3], ty[i+3], tz[i+3]},
                 xp={x[ofs], y[ofs], z[ofs]};
          /* Sub-tet volumes */
          double Va=dotProduct(crossProduct(minus(dp,bp), minus(cp,bp)),
                                            minus(xp,bp))*invV0[j];
          double Vb=dotProduct(crossProduct(minus(cp,ap), minus(dp,ap)),
                                            minus(xp,ap))*invV0[j];
          double Vc=dotProduct(crossProduct(minus(dp,ap), minus(bp,ap)),
                                            minus(xp,ap))*invV0[j];
          //double Vd=dotProduct(crossProduct(minus(bp,ap), minus(cp,ap)),minus(xp,ap))*invV0[j];
          double Vd=1.0-Va-Vb-Vc;
          /*If the point is inside the triangle */
          /*A negative threshold is necessary if the interpolation point is
            on a tetrahedron edge and due to the numerical errors buried
            in Va..Vd */
          if ((Va >= -1e-13) && (Vb >= -1e-13) && (Vc >= -1e-13) && (Vd >= -1e-13)){
            switch (elementType){
              case P0:
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuIm[nArg] != NULL){
                    *(wIm[nArg]+ofs)=*(tuIm[nArg]+j);
                    *(wRe[nArg]+ofs)=*(tuRe[nArg]+j);
                  }
                  else{
                    *(wRe[nArg]+ofs)=*(tuRe[nArg]+j);
                  }
                }
              break;
              case P1:
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuIm[nArg] != NULL){
                    *(wIm[nArg]+ofs)=P1(Va,Vb,Vc,Vd,
                                        *(tuIm[nArg]+i),*(tuIm[nArg]+i+1),
                                        *(tuIm[nArg]+i+2),*(tuIm[nArg]+i+3));
                    *(wRe[nArg]+ofs)=P1(Va,Vb,Vc,Vd,
                                        *(tuRe[nArg]+i),*(tuRe[nArg]+i+1),
                                        *(tuRe[nArg]+i+2),*(tuRe[nArg]+i+3));
                  }
                  else{
                    *(wRe[nArg]+ofs)=P1(Va,Vb,Vc,Vd,
                                        *(tuRe[nArg]+i),*(tuRe[nArg]+i+1),
                                        *(tuRe[nArg]+i+2),*(tuRe[nArg]+i+3));
                  }
                }
              break;
              case P2:
                {
                mwSize k=j*10;
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuIm[nArg] != NULL){
                    *(wIm[nArg]+ofs)=P2(Va,Vb,Vc,
                                        *(tuIm[nArg]+k),*(tuIm[nArg]+k+1),
                                        *(tuIm[nArg]+k+2),*(tuIm[nArg]+k+3),
                                        *(tuIm[nArg]+k+4),*(tuIm[nArg]+k+5),
                                        *(tuIm[nArg]+k+6),*(tuIm[nArg]+k+7),
                                        *(tuIm[nArg]+k+8),*(tuIm[nArg]+k+9));
                    *(wRe[nArg]+ofs)=P2(Va,Vb,Vc,
                                        *(tuRe[nArg]+k),*(tuRe[nArg]+k+1),
                                        *(tuRe[nArg]+k+2),*(tuRe[nArg]+k+3),
                                        *(tuRe[nArg]+k+4),*(tuRe[nArg]+k+5),
                                        *(tuRe[nArg]+k+6),*(tuRe[nArg]+k+7),
                                        *(tuRe[nArg]+k+8),*(tuRe[nArg]+k+9));
                  }
                  else{
                    *(wRe[nArg]+ofs)=P2(Va,Vb,Vc,
                                        *(tuRe[nArg]+k),*(tuRe[nArg]+k+1),
                                        *(tuRe[nArg]+k+2),*(tuRe[nArg]+k+3),
                                        *(tuRe[nArg]+k+4),*(tuRe[nArg]+k+5),
                                        *(tuRe[nArg]+k+6),*(tuRe[nArg]+k+7),
                                        *(tuRe[nArg]+k+8),*(tuRe[nArg]+k+9));
                  }
                }
                }
              break;
              default:
                   //poormans fall through
              break;
            }
            doSearchTet=false;
          }
        }
        i=i+4;
        j=j+1;
      }
    }
  }
  mxFree(invV0);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

  const mwSize nArgMax=100;

  feElement elementType=P1;
  double *x, *y, *z;
  mwSize nColX, mRowX, nColY, mRowY, nColZ, mRowZ;
  double *tx, *ty, *tz;
  mwSize nColTX, mRowTX, nColTY, mRowTY, nColTZ, mRowTZ;
  double *tuRe[nArgMax], *tuIm[nArgMax];
  mwSize nColTU[nArgMax], mRowTU[nArgMax];
  double *wRe[nArgMax], *wIm[nArgMax];

  if ((nrhs < 7) || (nrhs > nArgMax)){
     mexErrMsgTxt("Number of input arguments not plausible");
  }

  mwSize nDim=(mwSize)(nrhs-6);
  //mexPrintf("nnDim: %i\n",nDim);

  /*Meshgrid x,y */
  x=mxGetPr(prhs[0]);
  nColX=(mwSize)mxGetN(prhs[0]);
  mRowX=(mwSize)mxGetM(prhs[0]);
  y=mxGetPr(prhs[1]);
  nColY=(mwSize)mxGetN(prhs[1]);
  mRowY=(mwSize)mxGetM(prhs[1]);
  z=mxGetPr(prhs[2]);
  nColZ=(mwSize)mxGetN(prhs[2]);
  mRowZ=(mwSize)mxGetM(prhs[2]);
  if (!((nColX == nColY) && (mRowX == mRowY) && (nColX == nColZ) && (mRowX == mRowZ))){
      mexErrMsgTxt("Arguments 1,2,3 must have same dimensions");
  }
  mwSize nX=nColX;
  mwSize nY=mRowX;

  /*Triangle Mesh tx, ty */
  tx=mxGetPr(prhs[3]);
  nColTX=(mwSize)mxGetN(prhs[3]);
  mRowTX=(mwSize)mxGetM(prhs[3]);
  ty=mxGetPr(prhs[4]);
  nColTY=(mwSize)mxGetN(prhs[4]);
  mRowTY=(mwSize)mxGetM(prhs[4]);
  tz=mxGetPr(prhs[5]);
  nColTZ=(mwSize)mxGetN(prhs[5]);
  mRowTZ=(mwSize)mxGetM(prhs[5]);
  if (!((mRowTX == 4) && (mRowTY == 4) && (mRowTZ == 4))){
    mexErrMsgTxt("Arguments 4,5,6 must have 4 rows");
  }
  if (!((nColTX == nColTY) && (nColTX == nColTZ))){
    mexErrMsgTxt("Arguments 4,5,6 must have same number of columns");
  }
  mwSize nTet=nColTX;

  /*FE data tu */
  for (mwSize i=0; i<nDim; i++){
    nColTU[i]=(mwSize)mxGetN(prhs[i+6]);
    mRowTU[i]=(mwSize)mxGetM(prhs[i+6]);
    if (mxIsComplex(prhs[i+6])){
      tuRe[i]=mxGetPr(prhs[i+6]);
      tuIm[i]=mxGetPi(prhs[i+6]);
    }
    else{
      tuRe[i]=mxGetPr(prhs[i+6]);
      tuIm[i]=NULL;
    }
  }
  mwSize nDoF=mRowTU[0];

  for (mwSize i=0; i<nDim; i++){
    if (!((nDoF == mRowTU[i]) && (nTet == nColTU[i]))){
      mexErrMsgTxt("Data Arguments must be all equal and have nT datapoints");
    }
  }

  //mexPrintf("nXcol:%i nYrow:%i nT:%i nDoF:%i\n",nX,nY,nTet,nDoF);

  switch (nDoF){
    case 1:
      elementType=P0;
    break;
    case 4:
      elementType=P1;
    break;
    case 10:
      elementType=P2;
    break;
    default:
      mexErrMsgTxt("Could not determine Type of Lagrangian Finite Element");
    break;
  }

  if (nlhs != nDim){
    mexErrMsgTxt("Same dimension for output as well as input");
  }

  for (mwSize i=0; i<nDim; i++){
    if (tuIm[i] != NULL){
      plhs[i]=mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
      wRe[i]=mxGetPr(plhs[i]);
      wIm[i]=mxGetPi(plhs[i]);
    }
    else{
      plhs[i]=mxCreateDoubleMatrix(nY, nX, mxREAL);
      wRe[i]=mxGetPr(plhs[i]);
      wIm[i]=NULL;
    }
  }

  fftet2gridfast (x,y,z,tx,ty,tz,tuRe,tuIm,wRe,wIm,nDim,nTet,
                  nX,nY,elementType);
}
