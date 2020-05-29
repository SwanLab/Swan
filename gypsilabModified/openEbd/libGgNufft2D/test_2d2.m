nj = 10000;
ms = 210;
mt = 210;

xj = sort((rand(nj,1)*2-1)*pi);
yj = sort((rand(nj,1)*2-1)*pi);
fk = randn(ms,mt)+1i*randn(ms,mt);

eps = 1e-12;


iflag = +1;

tic
cj = dirft2d2(nj,xj,yj,iflag,ms,mt,fk);
toc

tic
cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk);
toc

abs_error=norm(cj-cj1,2)
rel_error=norm(cj-cj1,2)/norm(cj,2)


iflag = -1;

tic
cj = dirft2d2(nj,xj,yj,iflag,ms,mt,fk);
toc

tic
cj1 = nufft2d2(nj,xj,yj,iflag,eps,ms,mt,fk);
toc

abs_error=norm(cj-cj1,2)
rel_error=norm(cj-cj1,2)/norm(cj,2)
