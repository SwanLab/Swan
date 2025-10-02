nx=200; ny=150; nxy=nx*ny;
eps=10;
k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
kperp=[-ky;kx];
alpha=4;

Dx=sparse(nxy,nxy); Dy=sparse(nxy,nxy);
ep=1e-5;
X=ones(ny,1)*[1:nx]; X=reshape(X',1,nx*ny)';
Y=[1:ny]'*ones(1,nx); Y=reshape(Y',1,nx*ny)';
%fidi=0.5*[-1 0 1];
fidi=[-1 1 0];
Dx=spdiags(ones(nxy,1)*fidi,-1:1,nxy,nxy);
Dx([1:nx:nxy],:)=0; Dx([nx:nx:nxy],:)=0;
Dy=spdiags(ones(nxy,1)*fidi,[-nx 0 nx],nxy,nxy);
Dy([1:nx],:)=0; Dy([nxy-nx+1:nxy],:)=0;

% Initialization
%U=real(X+Y<0.35*(nx+ny))+real(X+Y>0.65*(nx+ny)); % band
%U=real((X-nx/2).^2+(Y-ny/2).^2<min(nx,ny)^2/16); % ball
U=real(abs(X-nx/2)<min(nx,ny)/6).*real(abs(Y-ny/2)<min(nx,ny)/6); % square
vol=sum(U);

% Main loop
disp('main loop');
%noDx=norm(full(Dx)); noDy=norm(full(Dy));
sigmacp=0.25; taucp=0.25; thetacp=1;
taucpe=taucp/eps^2;
for iopt=1:20
% Chambolle Pock
%V=U; Vbar=U;
V=zeros(nxy,1); Vbar=zeros(nxy,1);
Wx=zeros(nxy,1); Wy=zeros(nxy,1);
for kcp=1:1000
    xix=Wx+sigmacp*Dx*Vbar; xiy=Wy+sigmacp*Dy*Vbar;
    xi=[xix,xiy];
    W1=xi/(1+sigmacp);
    Wx1=W1(:,1); Wy1=W1(:,2);
    xik=xi(:,1)*k(1)+xi(:,2)*k(2);
    xikperp=xi(:,1)*kperp(1)+xi(:,2)*kperp(2);
    W2k=((alpha^4*(1+sigmacp)+sigmacp^2*(alpha^2-1))/(alpha^2*(alpha^2+sigmacp)*(1+sigmacp)))*xik-((sigmacp*sqrt(alpha^2-1))/(alpha^2*(1+sigmacp)))*abs(xikperp);
    sig=real(xikperp>0)-real(xikperp<0);
    W2kperp=((alpha^2+sigmacp)*xikperp-sigmacp*sqrt(alpha^2-1)*sig.*xik)/(alpha^2*(1+sigmacp));
    W2=W2k*k'+W2kperp*kperp';
    Wx2=W2(:,1); Wy2=W2(:,2);
    noxi=sqrt(xi(:,1).^2+xi(:,2).^2);
    case1=real(alpha*xik-noxi<0);
    case2=real(alpha*xik-noxi>0);
    Wx=case1.*Wx1+case2.*Wx2; Wy=case1.*Wy1+case2.*Wy2;
    Vold=V;
    zeta=Vold-taucp*Dx'*Wx-taucp*Dy'*Wy;
    V=(zeta+taucpe*U)/(1+taucpe);
    Vbar=V+thetacp*(V-Vold);
end
Ut=reshape(U,nx,ny)';
figure(1); clf; surf(Ut);
UC=0.8*flipud(-Ut+1);
C=UC; C(:,:,2)=UC; C(:,:,3)=UC;
figure(3); clf; image(C); set(gca,'xtick',[]); set(gca,'ytick',[]); axis image; 

Vt=reshape(V,nx,ny)';
figure(2); clf; surf(Vt);
%pause;
% Bissection
lambdamin=0; lambdamax=1;
for idic=1:100
    lambda=(lambdamin+lambdamax)/2;
    Utest=real(V>lambda);
    voltest=sum(Utest);
    if (voltest>vol)
        lambdamin=lambda;
    else
        lambdamax=lambda;
    end
end
U=Utest;
sum(U)
end