function dropletperopt3

nx=50; ny=50; 
%nxy=nx*ny;
k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
kperp=[-ky;kx];
alpha=4;

Lx = 1; Ly = 1;


m = QuadMesh(Lx,Ly,nx,ny);

%D = createDerivative(nx,ny,nxy,dx,dy);

% Initialization
fracx = 1/3;
fracy = 1/3;
cx = Lx/2; 
cy = Ly/2;                 
w  = fracx * Lx; h = fracy * Ly;

s.xSide = w;
s.ySide = h;
s.xCoorCenter = cx;
s.yCoorCenter = cy;
s.type = 'Rectangle';
g = GeometricalFunction(s);
ls = g.computeLevelSetFunction(m);

sU.backgroundMesh = m;
sU.boundaryMesh   = m.createBoundaryMesh;
uMesh             = UnfittedMesh(sU);
uMesh.compute(ls.fValues);
cFun = CharacteristicFunction.create(uMesh);
s.filterType = 'LUMP';
s.mesh  = m;
s.trial = createRho(m);
f = Filter.create(s);
d = f.compute(cFun,2);

rho0 = d;
vol= Integrator.compute(rho0,m,2);

%noDx=norm(full(Dx)); noDy=norm(full(Dy));

tauF = 0.025; 
tauG = 0.25; 
thetaRel = 1;
ep=1e-5;
eps=10; 
taucpe=tauG/eps^2;

%proxF = @(z)  proximalDroplet(z,tauF,k,alpha,ep);
proxF = @(z)  proximalEllipse(z,tauF,alpha,m);
proxG = @(rho,chi) proximalL2Projection(rho,chi,taucpe);
%grad = @(u) Grad(u);
%div  = @(z) Div(z);

rho = rho0;
z0 = createSigmaFunction(m);
for iopt=1:20
[rho,z] = PerimeterMinimization(rho,z0,proxF,proxG,tauF,tauG,thetaRel);
rho = ProjectToVolumeConstraint(rho,vol);
close all
rho.plot()
drawnow
volCon = computeVolumeConstraint(rho,vol)
end

end

function volCon = computeVolumeConstraint(rho,vol)
volCon = Integrator.compute(rho,rho.mesh,2)/vol - 1;
end


function s = proximalEllipse(z,tau,alpha,m)
A = [1 0; 0 1];
I = eye(2);
r = alpha^2/tau;
dm = createTensorFunction(A+r*I,m);
s = createSigmaFunction(m);
M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
F  = IntegrateRHS(@(v) DP(v,z),s,m,'Domain');
sV = full(M\F);
sV = reshape(sV,[s.ndimf,m.nnodes])';
s.setFValues(sV);
end

function s = createSigmaFunction(m)
s = LagrangianFunction.create(m,2,'P1');
end

function rho = createRho(m)
rho = LagrangianFunction.create(m,1,'P1');
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

function rho = proximalL2Projection(rho,Chi,tauG)
rho = project(rho + tauG*Chi,'P1');
end

function [u,z] = solveWithChambollePockAlgorithm(u0,z0,proxF,proxG,tauF,tauG,thetaRel)
u  = u0; 
uN = u0;
z   = z0; 
for kcp=1:10
    z      = proxF(z + tauF.*Grad(uN));
    uOld   = u;
    u      = proxG(u - tauG*Divergence(z));
    uN     = project(u + thetaRel.*(u - uOld),'P1');
end
end

function [u,z] = PerimeterMinimization(u0,z0,proxF,proxGX,tauF,tauG,thetaRel)
proxG = @(rhoX) proxGX(rhoX,u0);   
[u,z] = solveWithChambollePockAlgorithm(u0,z0,proxF,proxG,tauF,tauG,thetaRel);
end

function Af = createTensorFunction(A,m)
op = @(xV) evaluate(xV,A,m);
Af = DomainFunction.create(op,m,2);
end

function fV = evaluate(xV,A,m)
nGauss = size(xV,2);
nElem  = m.nelem;
fV = repmat(A,[1 1 nGauss nElem]);
end





function plot(rho0,nx,ny)
Ut=reshape(rho0,nx,ny)';
figure(1); clf; surf(Ut);
UC=0.8*flipud(-Ut+1);
C=UC; C(:,:,2)=UC; C(:,:,3)=UC;
figure(3); clf; image(C); set(gca,'xtick',[]); set(gca,'ytick',[]); axis image; 
end

function U = ProjectToVolumeConstraint(u,vol)
U = createRho(u.mesh);
lLb=0; lUb=1;
for idic=1:100
    lam = (lLb+lUb)/2;
    U.setFValues(u.fValues > lam);
    constraint = computeVolumeConstraint(U,vol);
    if (constraint > 0)
        lLb =lam;
    else
        lUb =lam;
    end
end
U.setFValues(real(U.fValues));
end