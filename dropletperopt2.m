function dropletperopt2

nx=200; ny=150; nxy=nx*ny;
k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
kperp=[-ky;kx];
alpha=4;


D = createDerivative(nx,nxy);

% Initialization
rho0 = inizalization(nx,ny);
vol=sum(rho0);

%noDx=norm(full(Dx)); noDy=norm(full(Dy));

tauF = 0.25; 
tauG = 0.25; 
thetaRel = 1;
ep=1e-5;
eps=10; 
taucpe=tauG/eps^2;

proxF = @(rho)      proximalDroplet(rho,tauF,k,alpha,ep);
proxG = @(rho,rho0) proximalL2Projection(rho,rho0,taucpe);
grad = @(u) Grad(D,u);
div  = @(z) Div(D,z);

rho = rho0;
z0  = zeros(nxy,2);
z = z0;
for iopt=1:20
[rho,z] = PerimeterMinimization(rho,z0,grad,div,proxF,proxG,tauF,tauG,thetaRel);
rho = ProjectToVolumeConstraint(rho,vol);
plotSurf(rho,nx,ny)
plot(rho,nx,ny)
volCon = computeVolumeConstraint(rho,vol)
end

end

function volCon = computeVolumeConstraint(rho,vol)
volCon = sum(rho)/vol - 1;
end

function U = inizalization(nx,ny)
[X,Y] = createSquareDomain(nx,ny);
%U=real(X+Y<0.35*(nx+ny))+real(X+Y>0.65*(nx+ny)); % band
%U=real((X-nx/2).^2+(Y-ny/2).^2<min(nx,ny)^2/16); % ball
U=real(abs(X-nx/2)<min(nx,ny)/6).*real(abs(Y-ny/2)<min(nx,ny)/6); % square
end

function [X,Y] = createSquareDomain(nx,ny)
X=ones(ny,1)*[1:nx];
X=reshape(X',1,nx*ny)';
Y=[1:ny]'*ones(1,nx); 
Y=reshape(Y',1,nx*ny)';
end

function D = createDerivative(nx,nxy)
    fidi = [-1 1 0];
    D = [spdiags(ones(nxy,1)*fidi,-1:1,nxy,nxy); ...
         spdiags(ones(nxy,1)*fidi,[-nx 0 nx],nxy,nxy)];
    D(1:nx,:) = 0; 
    D(nx:nx:end-nx,:) = 0; 
    D(end-nx+1:end,:) = 0;
end


function W = proximalDroplet(xi,sigmacp,k,alpha,ep)
W1 = xi/(1+sigmacp);
xik  = xi*k(:);
nxi2 = sum(xi.^2,2);
tA   = sigmacp/(alpha^2*(1+sigmacp));
tB   = sqrt((alpha^2-1)./(nxi2-xik.^2+ep));
deltaW = tA*(xi + (alpha^2-2)*xik.*k.' - tB.*((nxi2-2*xik.^2).*k.' + xik.*xi));
W      = W1 + (alpha*xik - sqrt(nxi2) > 0).*deltaW;
end

function rho = proximalL2Projection(rho,rho0,taucpe)
rho = (rho + taucpe*rho0)/(1+taucpe);
end

function [u,z] = solveWithChambollePockAlgorithm(u0,z0,Grad,Div,proxF,proxG,tauF,tauG,thetaRel)
u = u0; uN = u0;
z   = z0; 
for kcp=1:1000
    z      = proxF(z + tauF*Grad(uN));
    uOld   = u;
    u      = proxG(u - tauG*Div(z));
    uN     = u + thetaRel*(u - uOld);
end
end

function [u,z] = PerimeterMinimization(u0,z0,Grad,Div,proxF,proxGX,tauF,tauG,thetaRel)
proxG = @(rhoX) proxGX(rhoX,u0);   
[u,z] = solveWithChambollePockAlgorithm(u0,z0,Grad,Div,proxF,proxG,tauF,tauG,thetaRel);
end

function divZ = Div(D,z)
divZ = D.'*z(:);
end

function gradU = Grad(D,u)
nxy = size(u,1);
gradU = reshape(D*u,nxy,2);
end


function plotSurf(rho,nx,ny)
rhoT=reshape(rho,nx,ny)';
figure(2); clf; surf(rhoT);
end

function plot(rho0,nx,ny)
Ut=reshape(rho0,nx,ny)';
figure(1); clf; surf(Ut);
UC=0.8*flipud(-Ut+1);
C=UC; C(:,:,2)=UC; C(:,:,3)=UC;
figure(3); clf; image(C); set(gca,'xtick',[]); set(gca,'ytick',[]); axis image; 
end

function U = ProjectToVolumeConstraint(u,vol)
lLb=0; lUb=1;
for idic=1:100
    lam = (lLb+lUb)/2;
    uTest=  u > lam;
    constraint = computeVolumeConstraint(uTest,vol);
    if (constraint > 0)
        lLb =lam;
    else
        lUb =lam;
    end
end
U=real(uTest);
end