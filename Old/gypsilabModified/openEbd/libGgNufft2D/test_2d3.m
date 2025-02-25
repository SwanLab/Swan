nj = 10000;
nk = 21000;

xj = sort((rand(nj,1)*2-1)*pi);
yj = sort((rand(nj,1)*2-1)*pi);
cj = randn(nj,1)+1i*randn(nj,1);
sk = sort((rand(nk,1)*2-1)*pi);
tk = sort((rand(nk,1)*2-1)*pi);

eps=1e-6;


iflag = +1;

tic
fk = dirft2d3(nj,xj,yj,cj,iflag,nk,sk,tk);
toc

tic
fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

tic
fk = dirft2d3(nj,xj,yj,cj,iflag,nk,sk,tk);
toc

tic
fk1 = nufft2d3(nj,xj,yj,cj,iflag,eps,nk,sk,tk);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
