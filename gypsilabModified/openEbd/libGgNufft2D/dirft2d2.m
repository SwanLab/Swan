function cj=dirft2d2(nj,xj,yj,iflag,ms,mt,fk)
%DIRFT2D2: Direct (slow) computation of nonuniform FFT in R^2 - Type 2.
%
%  CJ = DIRFT2D2(NJ,XJ,YJ,IFLAG,MS,MT,FK);
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
%     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
%     mt     number of Fourier modes given  [ -mt/2: (mt-1)/2 ]
%     fk     Fourier coefficient values (complex *16 array)
%
%  Output parameters:
%
%     cj     output values (complex *16 array)
%
%

cj=zeros(nj,1)+1i*zeros(nj,1);

mex_id_ = 'dirft2d2(i int[x], i double[], i double[], io dcomplex[], i int[x], i int[x], i int[x], i dcomplex[])';
[cj] = nufft2d(mex_id_, nj, xj, yj, cj, iflag, ms, mt, fk, 1, 1, 1, 1);


