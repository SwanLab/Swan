/* Interpolates from a 3D tetrahedral mesh to a 2D grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2019-02-22
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
 *   IMPORTANT: This code is compatible with MATLAB version R2018 and later
 *   versions only (uses new the "Interleaved Complex API").
 *
 *   To build the source on Windows with (mingw) gcc execute:
 *
 *   setenv('MW_MINGW64_LOC','path to mingw');
 *   mex -R2018a fftet2gridfast.c
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
mxDouble max4(mxDouble i, mxDouble j, mxDouble k, mxDouble l){
  return max1(max1(i,j), max1(k,l));
}

static inline
mxDouble min4(mxDouble i, mxDouble j, mxDouble k, mxDouble l){
  return min1(min1(i,j), min1(k,l));
}

typedef struct{
  mxDouble i,j,k;
}Vector;

static inline
mxDouble dotProduct(Vector a, Vector b){
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
fftet2gridfast (mxDouble *x, mxDouble *y, mxDouble *z, mxDouble *tx, mxDouble *ty,
                mxDouble *tz, mxDouble **tuRe, mxComplexDouble **tuCplx, mxDouble **wRe,
                mxComplexDouble **wCplx, mwSize nOuts, mwSize nTet, mwSize nx,
                mwSize ny, feElement elementType){

  mxDouble *invV0=(mxDouble *)mxMalloc(nTet*sizeof(mxDouble));
  mxDouble init=mxGetNaN( );

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
        if (tuCplx[nArg] != NULL){
          (*(wCplx[nArg]+ofs)).real = init;
          (*(wCplx[nArg]+ofs)).imag = init;
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
          mxDouble Va=dotProduct(crossProduct(minus(dp,bp), minus(cp,bp)),
                                              minus(xp,bp))*invV0[j];
          mxDouble Vb=dotProduct(crossProduct(minus(cp,ap), minus(dp,ap)),
                                              minus(xp,ap))*invV0[j];
          mxDouble Vc=dotProduct(crossProduct(minus(dp,ap), minus(bp,ap)),
                                              minus(xp,ap))*invV0[j];
          //mxDouble Vd=dotProduct(crossProduct(minus(bp,ap), minus(cp,ap)),minus(xp,ap))*invV0[j];
          mxDouble Vd=1.0-Va-Vb-Vc;
          /*If the point is inside the triangle */
          /*A negative threshold is necessary if the interpolation point is
            on a tetrahedron edge and due to the numerical errors buried
            in Va..Vd */
          if ((Va >= -1e-13) && (Vb >= -1e-13) && (Vc >= -1e-13) && (Vd >= -1e-13)){
            switch (elementType){
              case P0:
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuCplx[nArg] != NULL){
                    (*(wCplx[nArg]+ofs)).imag=(*(tuCplx[nArg]+j)).imag;
                    (*(wCplx[nArg]+ofs)).real=(*(tuCplx[nArg]+j)).real;
                  }
                  else{
                    *(wRe[nArg]+ofs)=*(tuRe[nArg]+j);
                  }
                }
              break;
              case P1:
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuCplx[nArg] != NULL){
                    (*(wCplx[nArg]+ofs)).imag=P1(Va,Vb,Vc,Vd,
                                                 (*(tuCplx[nArg]+i)).imag,(*(tuCplx[nArg]+i+1)).imag,
                                                 (*(tuCplx[nArg]+i+2)).imag,(*(tuCplx[nArg]+i+3)).imag);
                    (*(wCplx[nArg]+ofs)).real=P1(Va,Vb,Vc,Vd,
                                                 (*(tuCplx[nArg]+i)).real,(*(tuCplx[nArg]+i+1)).real,
                                                 (*(tuCplx[nArg]+i+2)).real,(*(tuCplx[nArg]+i+3)).real);
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
                  if (tuCplx[nArg] != NULL){
                    (*(wCplx[nArg]+ofs)).imag=P2(Va,Vb,Vc,
                                                 (*(tuCplx[nArg]+k)).imag,(*(tuCplx[nArg]+k+1)).imag,
                                                 (*(tuCplx[nArg]+k+2)).imag,(*(tuCplx[nArg]+k+3)).imag,
                                                 (*(tuCplx[nArg]+k+4)).imag,(*(tuCplx[nArg]+k+5)).imag,
                                                 (*(tuCplx[nArg]+k+6)).imag,(*(tuCplx[nArg]+k+7)).imag,
                                                 (*(tuCplx[nArg]+k+8)).imag,(*(tuCplx[nArg]+k+9)).imag);
                    (*(wCplx[nArg]+ofs)).real=P2(Va,Vb,Vc,
                                                 (*(tuCplx[nArg]+k)).real,(*(tuCplx[nArg]+k+1)).real,
                                                 (*(tuCplx[nArg]+k+2)).real,(*(tuCplx[nArg]+k+3)).real,
                                                 (*(tuCplx[nArg]+k+4)).real,(*(tuCplx[nArg]+k+5)).real,
                                                 (*(tuCplx[nArg]+k+6)).real,(*(tuCplx[nArg]+k+7)).real,
                                                 (*(tuCplx[nArg]+k+8)).real,(*(tuCplx[nArg]+k+9)).real);
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
  mxDouble *x, *y, *z;
  mwSize nColX, mRowX, nColY, mRowY, nColZ, mRowZ;
  mxDouble *tx, *ty, *tz;
  mwSize nColTX, mRowTX, nColTY, mRowTY, nColTZ, mRowTZ;
  mxComplexDouble *tuCplx[nArgMax], *wCplx[nArgMax];
  mxDouble *tuRe[nArgMax], *wRe[nArgMax];
  mwSize nColTU[nArgMax], mRowTU[nArgMax];

  if ((nrhs < 7) || (nrhs > nArgMax)){
     mexErrMsgTxt("Number of input arguments not plausible");
  }

  mwSize nDim=(mwSize)(nrhs-6);
  mexPrintf("nDim: %i\n",nDim);

  /*Meshgrid x,y */
  x=mxGetDoubles(prhs[0]);
  nColX=(mwSize)mxGetN(prhs[0]);
  mRowX=(mwSize)mxGetM(prhs[0]);
  y=mxGetDoubles(prhs[1]);
  nColY=(mwSize)mxGetN(prhs[1]);
  mRowY=(mwSize)mxGetM(prhs[1]);
  z=mxGetDoubles(prhs[2]);
  nColZ=(mwSize)mxGetN(prhs[2]);
  mRowZ=(mwSize)mxGetM(prhs[2]);
  if (!((nColX == nColY) && (mRowX == mRowY) && (nColX == nColZ) && (mRowX == mRowZ))){
      mexErrMsgTxt("Arguments 1,2,3 must have same dimensions");
  }
  mwSize nX=nColX;
  mwSize nY=mRowX;

  /*Triangle Mesh tx, ty */
  tx=mxGetDoubles(prhs[3]);
  nColTX=(mwSize)mxGetN(prhs[3]);
  mRowTX=(mwSize)mxGetM(prhs[3]);
  ty=mxGetDoubles(prhs[4]);
  nColTY=(mwSize)mxGetN(prhs[4]);
  mRowTY=(mwSize)mxGetM(prhs[4]);
  tz=mxGetDoubles(prhs[5]);
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
      tuCplx[i] = mxGetComplexDoubles(prhs[i+6]);
      tuRe[i]=NULL;
    }
    else{
      tuRe[i] = mxGetDoubles(prhs[i+6]);
      tuCplx[i] = NULL;
    }
  }
  mwSize nDoF=mRowTU[0];

  for (mwSize i=0; i<nDim; i++){
    if (!((nDoF == mRowTU[i]) && (nTet == nColTU[i]))){
      mexErrMsgTxt("Data Arguments must be all equal and have nT datapoints");
    }
  }

  mexPrintf("nXcol:%i nYrow:%i nT:%i nDoF:%i\n",nX,nY,nTet,nDoF);

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
    if (tuCplx[i] != NULL){
      plhs[i]=mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
      wCplx[i]=mxGetComplexDoubles(plhs[i]);
      wRe[i]=NULL;
      mexPrintf("data %i is complex\n",i);
    }
    else{
      plhs[i]=mxCreateDoubleMatrix(nY, nX, mxREAL);
      wRe[i]=mxGetDoubles(plhs[i]);
      wCplx[i]=NULL;
      mexPrintf("data %i is real\n",i);
    }
  }

  fftet2gridfast (x,y,z,tx,ty,tz,tuRe,tuCplx,wRe,wCplx,nDim,nTet,
                  nX,nY,elementType);
}
