/* Interpolates from a 2D triangular mesh to a 2D cartesian or curved grid
 *
 * Author: Chloros2 <chloros2@gmx.de>
 * Created: 2019-02-22
 *
 *   This file is part of the ffmatlib which is hosted at
 *   https://github.com/samplemaker/freefem_matlab_octave_plot
 *
 *   [w1, ...] = fftri2gridfast (x, y, tx, ty, tu1, ...)
 *
 *   fftri2grid computes the function values w1, w2, ... over a mesh grid defined
 *   by the arguments x, y from a set of functions u1, u2, ... with values
 *   given on a triangular mesh tx, ty. The values are computed using first order
 *   or second order approximating basis functions (P1, P1b or P2 - Lagrangian Finite
 *   Elements). The function values w1, w2, ... are real if tu1, tu2, ... are real
 *   or complex if tu1, tu2, ... are complex. The mesh grid x, y can be cartesian
 *   or curved. fftri2grid returns NaNs when an interpolation point is outside the
 *   triangular mesh. fftri2gridfast.c is a runtime optimized mex implementation
 *   of the function fftri2grid.m.
 *
 *   IMPORTANT: This code is compatible with MATLAB version R2018 and later
 *   versions only (uses new the "Interleaved Complex API").
 *
 *   To build the source on Windows with (mingw) gcc execute:
 *
 *   setenv('MW_MINGW64_LOC','path to mingw');
 *   mex -R2018a fftri2gridfast.c
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

#include "mex.h"

#define min1(a,b) ( ((a)<(b)) ? (a):(b) )
#define max1(a,b) ( ((a)>(b)) ? (a):(b) )
#define min3(a,b,c) ( ((a)<(b)) ? (((a)<(c)) ? (a):(c)) : (((b)<(c)) ? (b):(c)) )
#define max3(a,b,c) ( ((a)>(b)) ? (((a)>(c)) ? (a):(c)) : (((b)>(c)) ? (b):(c)) )

#define P1(Aa,Ab,Ac,u1,u2,u3) (u1*Aa+u2*Ab+u3*Ac)

#define P1b(Aa,Ab,Ac,u1,u2,u3,u4) (u1*Aa+                              \
                                   u2*Ab+                              \
                                   u3*Ac+                              \
                                   9.0*(3.0*u4-(u1+u2+u3))*Aa*Ab*Ac)

#define P2(Aa,Ab,u1,u2,u3,u4,u5,u6) u1*Aa*(2.0*Aa-1.0)+               \
                                    u2*Ab*(2.0*Ab-1.0)+               \
                                    u3*(1.0-Aa-Ab)*(1.0-2.0*(Aa+Ab))+ \
                                    u4*4.0*Ab*(1.0-Aa-Ab)+            \
                                    u5*4.0*Aa*(1.0-Aa-Ab)+            \
                                    u6*4.0*Aa*Ab

typedef enum {
  P0,
  P1,
  P1b,
  P2
} feElement;

void
fftri2gridfast(mxDouble *x, mxDouble *y, mxDouble *tx, mxDouble *ty,
               mxDouble **tuRe, mxComplexDouble **tuCplx, mxDouble **wRe, mxComplexDouble **wCplx,
               mwSize nOuts, mwSize nTri, mwSize nx,  mwSize ny,
               feElement elementType){

  mxDouble *invA0=(mxDouble *)mxMalloc(nTri*sizeof(mxDouble));
  mxDouble init=mxGetNaN( );

  /*Triangle areas */
  mwSize j=0;
  mwSize nNodes=3*nTri;
  for (mwSize i=0; i<nTri; i++){
    invA0[i]=1.0/((ty[j+1]-ty[j+2])*(tx[j]-tx[j+2])+(tx[j+2]-tx[j+1])*(ty[j]-ty[j+2]));
    j=j+3;
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
      bool doSearchTri=true;
      while (doSearchTri && (i<nNodes)){
        /*If the point (X,Y) is outside a square defined by the triangle
          vertices, the point can not be within the triangle */
        bool presel=((x[ofs] <= max3(tx[i],tx[i+1],tx[i+2])) &&
                     (x[ofs] >= min3(tx[i],tx[i+1],tx[i+2])) &&
                     (y[ofs] <= max3(ty[i],ty[i+1],ty[i+2])) &&
                     (y[ofs] >= min3(ty[i],ty[i+1],ty[i+2])));
        /*Potential candiate - calculate Barycentric Coordinates */
        if (presel) {
          /*Sub-triangle areas */
          mxDouble Aa=((ty[i+1]-ty[i+2])*(x[ofs]-tx[i+2])+
                       (tx[i+2]-tx[i+1])*(y[ofs]-ty[i+2]))*invA0[j];
          mxDouble Ab=((ty[i+2]-ty[i])*(x[ofs]-tx[i+2])+
                       (tx[i]-tx[i+2])*(y[ofs]-ty[i+2]))*invA0[j];
          mxDouble Ac=1.0-Aa-Ab;
          /*If the point is inside the triangle */
          /*A negative threshold is necessary if the interpolation point is
            on a triangle edge and due to the numerical errors buried
            in Aa,b,c */
          if ((Aa >= -1e-13) && (Ab >= -1e-13) && (Ac >= -1e-13)){
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
                    (*(wCplx[nArg]+ofs)).imag=P1(Aa,Ab,Ac,
                                                 (*(tuCplx[nArg]+i)).imag,(*(tuCplx[nArg]+i+1)).imag,
                                                 (*(tuCplx[nArg]+i+2)).imag);
                    (*(wCplx[nArg]+ofs)).real=P1(Aa,Ab,Ac,
                                                 (*(tuCplx[nArg]+i)).real,(*(tuCplx[nArg]+i+1)).real,
                                                 (*(tuCplx[nArg]+i+2)).real);
                  }
                  else{
                    *(wRe[nArg]+ofs)=P1(Aa,Ab,Ac,
                                        *(tuRe[nArg]+i),*(tuRe[nArg]+i+1),
                                        *(tuRe[nArg]+i+2));
                  }
                }
              break;
              case P1b:
                {
                mwSize k=4*j;
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuCplx[nArg] != NULL){
                    (*(wCplx[nArg]+ofs)).imag=P1b(Aa,Ab,Ac,
                                                 (*(tuCplx[nArg]+k)).imag,(*(tuCplx[nArg]+k+1)).imag,
                                                 (*(tuCplx[nArg]+k+2)).imag,(*(tuCplx[nArg]+k+3)).imag);
                    (*(wCplx[nArg]+ofs)).real=P1b(Aa,Ab,Ac,
                                                 (*(tuCplx[nArg]+k)).real,(*(tuCplx[nArg]+k+1)).real,
                                                 (*(tuCplx[nArg]+k+2)).real,(*(tuCplx[nArg]+k+3)).real);
                  }
                  else{
                    *(wRe[nArg]+ofs)=P1b(Aa,Ab,Ac,
                                        *(tuRe[nArg]+k),*(tuRe[nArg]+k+1),
                                        *(tuRe[nArg]+k+2),*(tuRe[nArg]+k+3));
                  }
                }
              }
              break;
              case P2:
                {
                mwSize k=6*j;
                for (mwSize nArg=0; nArg<nOuts; nArg++){
                  if (tuCplx[nArg] != NULL){
                    (*(wCplx[nArg]+ofs)).imag=P2(Aa,Ab,
                                                 (*(tuCplx[nArg]+k)).imag,(*(tuCplx[nArg]+k+1)).imag,
                                                 (*(tuCplx[nArg]+k+2)).imag,(*(tuCplx[nArg]+k+3)).imag,
                                                 (*(tuCplx[nArg]+k+4)).imag,(*(tuCplx[nArg]+k+5)).imag);
                    (*(wCplx[nArg]+ofs)).real=P2(Aa,Ab,
                                                 (*(tuCplx[nArg]+k)).real,(*(tuCplx[nArg]+k+1)).real,
                                                 (*(tuCplx[nArg]+k+2)).real,(*(tuCplx[nArg]+k+3)).real,
                                                 (*(tuCplx[nArg]+k+4)).real,(*(tuCplx[nArg]+k+5)).real);
                  }
                  else{
                    *(wRe[nArg]+ofs)=P2(Aa,Ab,
                                        *(tuRe[nArg]+k),*(tuRe[nArg]+k+1),
                                        *(tuRe[nArg]+k+2),*(tuRe[nArg]+k+3),
                                        *(tuRe[nArg]+k+4),*(tuRe[nArg]+k+5));
                  }
                }
              }
              break;
              default:
                //poormans fall through
              break;
            }
            doSearchTri=false;
          }
        }
        i=i+3;
        j=j+1;
      }
    }
  }
  mxFree(invA0);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

  const mwSize nArgMax=100;

  feElement elementType=P1;
  mxDouble *x, *y;
  mwSize nColX, mRowX, nColY, mRowY;
  mxDouble *tx, *ty;
  mwSize nColTX, mRowTX, nColTY, mRowTY;
  mxComplexDouble *tuCplx[nArgMax], *wCplx[nArgMax];
  mxDouble *tuRe[nArgMax],*wRe[nArgMax];
  mwSize nColTU[nArgMax], mRowTU[nArgMax];

  if ((nrhs < 5) || (nrhs > nArgMax)){
    mexErrMsgTxt("Number of input arguments not plausible");
  }

  mwSize nDim=(mwSize)(nrhs-4);
  mexPrintf("nDim: %i\n",nDim);

  /*Meshgrid x,y */
  x=mxGetDoubles(prhs[0]);
  nColX=(mwSize)mxGetN(prhs[0]);
  mRowX=(mwSize)mxGetM(prhs[0]);
  y=mxGetDoubles(prhs[1]);
  nColY=(mwSize)mxGetN(prhs[1]);
  mRowY=(mwSize)mxGetM(prhs[1]);
  if (!((nColX == nColY) && (mRowX == mRowY))){
    mexErrMsgTxt("Arguments 1,2 must have same dimensions");
  }
  mwSize nX=nColX;
  mwSize nY=mRowX;

  /*Triangle Mesh tx, ty */
  tx=mxGetDoubles(prhs[2]);
  nColTX=(mwSize)mxGetN(prhs[2]);
  mRowTX=(mwSize)mxGetM(prhs[2]);
  ty=mxGetDoubles(prhs[3]);
  nColTY=(mwSize)mxGetN(prhs[3]);
  mRowTY=(mwSize)mxGetM(prhs[3]);
  if (!((mRowTX == 3) && (mRowTY == 3))){
    mexErrMsgTxt("Arguments 3,4 must have 3 rows");
  }
  if (!((nColTX == nColTY))){
    mexErrMsgTxt("Arguments 3,4 must have same number of columns");
  }
  mwSize nTri=nColTX;

  /*FE data tu */
  for (mwSize i=0; i<nDim; i++){
    nColTU[i]=(mwSize)mxGetN(prhs[i+4]);
    mRowTU[i]=(mwSize)mxGetM(prhs[i+4]);
    if (mxIsComplex(prhs[i+4])){
      tuCplx[i] = mxGetComplexDoubles(prhs[i+4]);
      tuRe[i]=NULL;
    }
    else{
      tuRe[i] = mxGetDoubles(prhs[i+4]);
      tuCplx[i] = NULL;
    }
  }
  mwSize nDoF=mRowTU[0];

  for (mwSize i=0; i<nDim; i++){
     if (!((nDoF == mRowTU[i]) && (nTri == nColTU[i]))){
       mexErrMsgTxt("Data Arguments must be all equal and have nT datapoints");
     }
  }

  mexPrintf("nXcol:%i nYrow:%i nT:%i nDoF:%i\n",nX,nY,nTri,nDoF);

  switch (nDoF){
    case 1:
      elementType=P0;
    break;
    case 3:
      elementType=P1;
    break;
    case 4:
      elementType=P1b;
    break;
    case 6:
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

  fftri2gridfast (x,y,tx,ty,tuRe,tuCplx,wRe,wCplx,nDim,nTri,nX,nY,
                  elementType);
}
