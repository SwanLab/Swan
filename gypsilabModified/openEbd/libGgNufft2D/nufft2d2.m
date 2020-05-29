function [cj,ier]=nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk)
%NUFFT2D2: Nonuniform FFT in R^2 - Type 2.
%
%  [CJ,IER] = NUFFT2D2(NJ,XJ,YJ,IFLAG,EPS,MS,MT,FK);
%
%     cj(j) = SUM   fk(k1,k2) exp(+/-i k1 xj(j)) exp(+/-i k2 yj(j)) 
%             k1,k2  
%                            for j = 1,...,nj
%
%     where -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2
%
%     If (iflag .ge.0) the + sign is used in the exponential.
%     If (iflag .lt.0) the - sign is used in the exponential.
%
%  Input parameters:
%
%     nj     number of output values   (integer)
%     xj,yj  location of output values (real *8 array)
%     iflag  determines sign of FFT (see above)
%     eps    precision request  (between 1.0e-15 and 1.0e-1)
%     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
%     mt     number of Fourier modes given  [ -mt/2: (mt-1)/2 ]
%     fk     Fourier coefficient values (complex *16 array)
%
%  Output parameters:
%
%     cj     output values (complex *16 array)
%     ier    error return code   
%            ier = 0  => normal execution.
%            ier = 1  => precision eps requested is out of range.
%
%

cj=zeros(nj,1)+1i*zeros(nj,1);
ier=0;

mex_id_ = 'nufft2d2f90(i int[x], i double[], i double[], io dcomplex[], i int[x], i double[x], i int[x], i int[x], i dcomplex[], io int[x])';
[cj, ier] = nufft2d(mex_id_, nj, xj, yj, cj, iflag, eps, ms, mt, fk, ier, 1, 1, 1, 1, 1, 1);


