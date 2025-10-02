function dropletperopt2

nx=200; ny=150; nxy=nx*ny;
k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
kperp=[-ky;kx];
alpha=4;


[Dx,Dy] = createDerivative(nx,nxy);

% Initialization
U = inizalization(nx,ny);
vol=sum(U);

% Main loop
disp('main loop');
%noDx=norm(full(Dx)); noDy=norm(full(Dy));

for iopt=1:20
V = PerimeterMinimization(U,nx,ny,nxy,alpha,k,Dx,Dy);
%pause;
U = ProjectToVolumeConstraint(V,vol);
sum(U)
end

end

function U = inizalization(nx,ny)
[X,Y] = createSquareDomain(nx,ny);
%U=real(X+Y<0.35*(nx+ny))+real(X+Y>0.65*(nx+ny)); % band
%U=real((X-nx/2).^2+(Y-ny/2).^2<min(nx,ny)^2/16); % ball
U=real(abs(X-nx/2)<min(nx,ny)/6).*real(abs(Y-ny/2)<min(nx,ny)/6); % square
end

function [X,Y] = createSquareDomain(nx,ny)
X=ones(ny,1)*[1:nx]; X=reshape(X',1,nx*ny)';
Y=[1:ny]'*ones(1,nx); Y=reshape(Y',1,nx*ny)';
end


function [Dx,Dy] = createDerivative(nx,nxy)
Dx=sparse(nxy,nxy); Dy=sparse(nxy,nxy);
%fidi=0.5*[-1 0 1];
fidi=[-1 1 0];
Dx=spdiags(ones(nxy,1)*fidi,-1:1,nxy,nxy);
Dx([1:nx:nxy],:)=0; Dx([nx:nx:nxy],:)=0;
Dy=spdiags(ones(nxy,1)*fidi,[-nx 0 nx],nxy,nxy);
Dy([1:nx],:)=0; Dy([nxy-nx+1:nxy],:)=0;
end










function V = PerimeterMinimization(U,nx,ny,nxy,alpha,k,Dx,Dy)
sigmacp=0.25; taucp=0.25; thetacp=1;ep=1e-5;eps=10;
taucpe=taucp/eps^2;
% Chambolle Pock
V=U; Vbar=U;
%V=zeros(nxy,1); Vbar=zeros(nxy,1);
Wx=zeros(nxy,1); Wy=zeros(nxy,1);
for kcp=1:1000
    xix=Wx+sigmacp*Dx*Vbar; xiy=Wy+sigmacp*Dy*Vbar;
    xi=[xix,xiy];
    W1=xi/(1+sigmacp);
    Wx1=W1(:,1); Wy1=W1(:,2);
    xik=xi(:,1)*k(1)+xi(:,2)*k(2);
    nxi2=xix.^2+xiy.^2;
    deltaWx2=(sigmacp/(alpha^2*(1+sigmacp)))*(xix+(alpha^2-2)*k(1)*xik-sqrt((alpha^2-1)./(nxi2-xik.^2+ep)).*((nxi2-2*xik.^2)*k(1)+xik.*xix));
    deltaWy2=(sigmacp/(alpha^2*(1+sigmacp)))*(xiy+(alpha^2-2)*k(2)*xik-sqrt((alpha^2-1)./(nxi2-xik.^2+ep)).*((nxi2-2*xik.^2)*k(2)+xik.*xiy));
    noxi=sqrt(nxi2);
    case2=real(alpha*xik-noxi>0);
    Wx=Wx1+case2.*deltaWx2; Wy=Wy1+case2.*deltaWy2;
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
end

function U = ProjectToVolumeConstraint(V,vol)
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
end