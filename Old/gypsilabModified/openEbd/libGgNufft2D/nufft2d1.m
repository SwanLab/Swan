function [fk,ier]=nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt)
%NUFFT2D1: Nonuniform FFT in R^2 - Type 1.
%
%  [FK,IER] = NUFFT2D1(NJ,XJ,YJ,CJ,IFLAG,EPS,MS,MT);
%
%                  1  nj
%     fk(k1,k2) = -- SUM cj(j) exp(+/-i k1 xj(j)) exp(+/-i k2 yj(j)) 
%                 nj j=1
%
%     for -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2
%
%     If (iflag .ge.0) the + sign is used in the exponential.
%     If (iflag .lt.0) the - sign is used in the exponential.
%
%  Input parameters:
%
%     nj     number of sources   (integer)
%     xj,yj  location of sources (real *8)
%
%            on interval [-pi,pi].
%
%     cj     strengths of sources (complex *16)
%     iflag  determines sign of FFT (see above)
%     eps    precision request  (between 1.0e-15 and 1.0e-1)
%     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
%     mt     number of Fourier modes computed (-mt/2 to (mt-1)/2 )
%                 
%  Output parameters:
%
%     fk     Fourier transform values (complex *16)
%     ier    error return code   
%            ier = 0  => normal execution.
%            ier = 1  => precision eps requested is out of range.
%
%

fk=zeros(ms,mt)+1i*zeros(ms,mt);
ier=0;

mex_id_ = 'nufft2d1f90(i int[x], i double[], i double[], i dcomplex[], i int[x], i double[x], i int[x], i int[x], io dcomplex[], io int[x])';
[fk, ier] = nufft2d(mex_id_, nj, xj, yj, cj, iflag, eps, ms, mt, fk, ier, 1, 1, 1, 1, 1, 1);


