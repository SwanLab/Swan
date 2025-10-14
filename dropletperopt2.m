function dropletperopt2

nx=300; ny=300; nxy=nx*ny;
k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
kperp=[-ky;kx];
alpha=4;

Lx = 1; Ly = 1;
dx = Lx/(nx-1); 
dy = Ly/(ny-1);

D = createDerivative(nx,ny,nxy,dx,dy);

% Initialization
chi0 = inizalization(nx,ny,Lx,Ly);
vol  =sum(chi0);

%noDx=norm(full(Dx)); noDy=norm(full(Dy));

tauF = 0.25; 
tauG = 0.25; 
thetaRel = 1;
ep=1e-5;
eps=10; 
taucpe=tauG/eps^2;

%proxF = @(z)  proximalDroplet(z,tauF,k,alpha,ep);
proxF = @(z)  proximalEllipse(z,tauF,alpha);
f     = 
df
g
dg

proxG = @(rho,rho0) proximalL2Projection(rho,rho0,taucpe);
grad = @(u) Grad(D,u);
div  = @(z) Div(D,z);

chi = chi0;
rho = chi0;
z0  = zeros(nxy,2);
z   = z0;
for iopt=1:20
[rho,z] = PerimeterComputation(chi,rho,z0,grad,div,proxF,proxG,tauF,tauG,thetaRel);
chi = ProjectToVolumeConstraint(rho,vol);
plotSurf(chi,nx,ny)
plot(chi,nx,ny)
volCon = computeVolumeConstraint(chi,vol)
end

end

function volCon = computeVolumeConstraint(chi,vol)
volCon = sum(chi)/vol - 1;
end

function U = inizalization(nx,ny,Lx,Ly)
[X,Y] = createSquareDomain(nx,ny,Lx,Ly);
%U=real(X+Y<0.35*(nx+ny))+real(X+Y>0.65*(nx+ny)); % band
%U=real((X-nx/2).^2+(Y-ny/2).^2<min(nx,ny)^2/16); % ball
%U=real(abs(X-nx/2)<min(nx,ny)/6).*real(abs(Y-ny/2)<min(nx,ny)/6); % square
    fracx = 1/3;
    fracy = 1/3;
    cx = Lx/2; cy = Ly/2;                 % center
    w  = fracx * Lx; h = fracy * Ly;      % sizes
    U  = double(abs(X-cx)<=w/2 & abs(Y-cy)<=h/2);

end

function [X,Y] = createSquareDomain(nx,ny,Lx,Ly)
[XX,YY] = meshgrid(linspace(0,Lx,nx), linspace(0,Ly,ny));
X = XX(:);
Y = YY(:);
end
 
function D = createDerivative(nx,ny,nxy,Lx,Ly)
    dx=1;Lx/(nx-1); dy=1;Ly/(ny-1);
    fidi = [-1 1 0];
    D = [spdiags(ones(nxy,1)*fidi*dx,-1:1,nxy,nxy); ...
         spdiags(ones(nxy,1)*fidi*dy,[-nx 0 nx],nxy,nxy)];
    D(1:nx,:) = 0; 
    D(nx:nx:end-nx,:) = 0; 
    D(end-nx+1:end,:) = 0;
end


function [s,J] = proximalEllipse(z,tau,alpha)
A = [1 0; 0 1];
I = eye(2);
r = alpha^2/tau;
invM = inv((A+r*I));
s = r*z*invM.';
J = 0.5*sum(sum((z*invM.').*z));
end


function s = proximalDroplet(z,tau,k,alpha,ep)
s1 = z/(1+tau);
zk  = z*k(:);
z2 = sum(z.^2,2);
tA   = tau/(alpha^2*(1+tau));
tB   = sqrt((alpha^2-1)./(z2-zk.^2+ep));
deltaS = tA*(z + (alpha^2-2)*zk.*k.' - tB.*((z2-2*zk.^2).*k.' + zk.*z));
s      = s1 + (alpha*zk - sqrt(z2) > 0).*deltaS;
end

function rhoN = proximalL2Projection(rho,Chi,tauG)
rhoN = (rho + tauG*Chi)/(1+tauG);
end

function J = L2ProjectionPrimal(rho,Chi)
J = 0.5*(rho-Chi)'*(rho-Chi);
end

function J = L2ProjectionDual(rho,Chi)
J = 0.5*(rho-Chi)'*(rho-Chi);
end

function J = EllipseGradientPrimal(txi)
A = [1 0; 0 1];
J = 0.5*sum(sum((txi*A.').*txi));
end


function J = EllipseGradientDual(txi)
A = [1 0; 0 1];
invA = inv(A);
J = 0.5*sum(sum((txi*invA.').*txi));
end

function [u,z] = solveWithChambollePockAlgorithm(u0,z0,Grad,Div,proxF,proxG,tauF,tauG,thetaRel)
u  = u0; 
uN = u0;
z   = z0; 
for kcp=1:2000
    [z,fd] = proxF(z + tauF*Grad(uN));
    uOld   = u;
    [u,gd] = proxG(u - tauG*Div(z));
    uN     = u + thetaRel*(u - uOld);
end
end

function [u,z] = PerimeterComputation(chi,u0,z0,Grad,Div,proxF,proxGX,tauF,tauG,thetaRel)
proxG = @(rho) proxGX(rho,chi);   
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

function chi = ProjectToVolumeConstraint(rho,vol)
lLb=0; lUb=1;
for idic=1:100
    lam = (lLb+lUb)/2;
    uTest=  rho > lam;
    constraint = computeVolumeConstraint(uTest,vol);
    if (constraint > 0)
        lLb =lam;
    else
        lUb =lam;
    end
end
chi=real(uTest);
end