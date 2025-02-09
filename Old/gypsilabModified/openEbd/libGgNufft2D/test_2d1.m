nj = 10000;
ms = 210;
mt = 210;

xj = sort((rand(nj,1)*2-1)*pi);
yj = sort((rand(nj,1)*2-1)*pi);
cj = randn(nj,1)+1i*randn(nj,1);

eps=1e-12;


iflag = +1;

tic
fk = dirft2d1(nj,xj,yj,cj,iflag,ms,mt);
toc

tic
fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

tic
fk = dirft2d1(nj,xj,yj,cj,iflag,ms,mt);
toc

tic
fk1 = nufft2d1(nj,xj,yj,cj,iflag,eps,ms,mt);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
