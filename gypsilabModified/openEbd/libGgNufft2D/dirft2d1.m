function fk=dirft2d1(nj,xj,yj,cj,iflag,ms,mt)
%DIRFT2D1: Direct (slow) computation of nonuniform FFT in R^2 - Type 1.
%
%  FK = DIRFT2D1(NJ,XJ,YJ,CJ,IFLAG,MS,MT);
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
%     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
%     mt     number of Fourier modes computed (-mt/2 to (mt-1)/2 )
%                 
%  Output parameters:
%
%     fk     Fourier transform values (complex *16)
%
%

fk=zeros(ms,mt)+1i*zeros(ms,mt);

mex_id_ = 'dirft2d1(i int[x], i double[], i double[], i dcomplex[], i int[x], i int[x], i int[x], io dcomplex[])';
[fk] = nufft2d(mex_id_, nj, xj, yj, cj, iflag, ms, mt, fk, 1, 1, 1, 1);


