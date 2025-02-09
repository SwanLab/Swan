%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal, Marc Bakry (c) 2017-2019.                 |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             marc.bakry@polytechnique.edu                               |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtHmxVibroSlab.m                             |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
tol = 1e-3
X0  = [0 0 1]
e   = 0.5;
f   = 900;

% Exterior domain (water)
rho0 = 1000;         % density      (kg.m3)
c0   = 1500;         % celerity     (m.s-1)
k0   = 2*pi/c0*f;    % wave-number  (m-1)
lam0 = c0/f;         % wave-length  (m)
    
% Interior domain (different from water)
rhoS = 2*rho0;                                 % density                             (kg.m3)
cL   = 2*c0;                                   % celerity of longitudinal waves      (m.s-1)
kL   = 2*pi*f/cL;                              % wave-number of longitudinal waves   (m-1)
lamL = real(cL)/f;                             % wavelength of longitudinal waves    (m)
cT   = 0;                                      % celerity of transverse waves        (m.s-1)
kT   = 2*pi.*f/cT;                             % wave-number of transverse waves     (m-1)
lamT = real(cT)/f;                             % wavelength of transverse waves 

% Minimum wavelength
tmp  = [lam0,lamL,lamT];
lmin = min(tmp(tmp>0));

% Slab mesh
L    = 20 * lmin;            % 20 wavelength to simulate infinite slab
nx   = ceil(L/lmin * 6)+1;   % 6 node per wavelength for L
ny   = ceil(L/lmin * 6)+1;   % 6 node per wavelength for L
nz   = ceil(e/lmin * 12)+1;  % 12 node per wavelength for e
N    = nx * ny * nz;         % Total number of nodes
mesh = mshCube(N,[L L e])

% Boundary
bound = mesh.bnd

% Admissible frequency
stp  = mesh.stp;
kmax = 1/stp(2);
if (k0 > kmax) || (kL > kmax/2) || ((kT>kmax/2) && ~isinf(kT))
%     error('frequency is too high for mesh discretization')
end

% Radiative mesh (fixed number of nodes)
radiatx     = mshSquare(1e3,[L L]);
radiaty     = radiatx;
radiatx.vtx = radiatx.vtx(:,[1 3 2]);
radiaty.vtx = radiaty.vtx(:,[3 1 2]);
radiat      = union(radiatx,radiaty);

% Measurement points for trans and refl coeff (1 wavelenth from bound)
Xmes = [0 0 -e/2-lmin ; 0 0 e/2+lmin];

% Cut-off function (50% full, 10% decrease)
cutoffx = vibsCutoff(1,L/5,L/10);
cutoffy = vibsCutoff(2,L/5,L/10);

% Green kernel function
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k0);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k0);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k0);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k0);
G0          = '[1/r]';
gradyG0     = 'grady[1/r]';
cteGxy      =  1/(4*pi);
cteG0       =  1/(4*pi);
    
% Plane wave function
PW         = @(X) exp(1i*k0*X*X0') .* cutoffx(X) .* cutoffy(X);
gradxPW{1} = @(X) 1i*k0*X0(1) .* PW(X);
gradxPW{2} = @(X) 1i*k0*X0(2) .* PW(X);
gradxPW{3} = @(X) 1i*k0*X0(3) .* PW(X);

% Coupling coeff for Brackage-Werner simulation
beta = 1i*k0;

% Graphical representation
figure
plot(bound,'r')
hold on
plot(radiat,log(abs(PW(radiat.vtx))))
% plotNrm(bound,'w')
plot(msh(Xmes),'k')
axis equal
colorbar
view(45,10)

% Quadrature and finite elements (volumn)
omega = dom(mesh,4);
U     = fem(mesh,'P1');

% Quadrature and finite elements (boundary)
sigma = dom(bound,3);
u     = fem(bound,'P1');

% Left-hand side
tic
[A,B,C,D] = vibsHmxBlockOperator(omega,U,sigma,u,cL,cT,rhoS,c0,rho0,f,tol);
toc

% Add dirichlet condition to x and y unknows (penalization)
A(sub2ind(size(A),1:2*length(U),1:2*length(U))) = 1e15;

% Right-hand side
V    = cell(4,1);
V{1} = - integral(sigma,ntimes(U,1),PW);
V{2} = - integral(sigma,ntimes(U,2),PW);
V{3} = - integral(sigma,ntimes(U,3),PW);
V{4} = integral(sigma,ntimes(u),gradxPW);

% Resolution with Schur complement
Fa     = decomposition(A);
LHS    = @(V) D*V - C*(Fa \ (B.Ml*(B.Mr*V)) );
RHS    = V{end} - C*(Fa \ cell2mat(V(1:end-1)) );
mu     = mgcr(LHS,RHS,[],tol,100);
lambda = beta*mu;

% Measure of refexive and transmitted coeff
tic
Pmes = cteGxy .* integral(Xmes,sigma,Gxy,u)*lambda - ...
    cteGxy .* integral(Xmes,sigma,gradyGxy,ntimes(u))*mu;
Pmes(2) = Pmes(2) + PW(Xmes(2,:));
toc

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy
tic
Srad = cteGxy .* integral(radiat.vtx,sigma,Gxy,u,tol);
Srad = Srad + cteG0 .* regularize(radiat.vtx,sigma,G0,u);

% Finite element radiative operator --> \int_Sy dnyG(x,y) psi(y) dy
Drad = cteGxy .* integral(radiat.vtx,sigma,gradyGxy,ntimes(u),tol);
Drad = Drad + cteG0 .* regularize(radiat.vtx,sigma,gradyG0,ntimes(u));
toc

% Domain solution
Psca = Srad*lambda - Drad*mu;
Pinc = PW(radiat.vtx);
Ptot = Psca + Pinc;

% Graphical repesentation
figure
subplot(1,2,1)
plot(bound)
hold on
plot(radiat,abs(Psca))
plot(msh(Xmes),abs(Pmes))
title('Scattered')
axis equal
colorbar
view(45,10)

subplot(1,2,2)
plot(bound)
hold on
plot(radiat,abs(Ptot))
plot(msh(Xmes),abs(Pmes))
title('Total')
axis equal
colorbar
view(45,10)

% Analytical solution
tic
c     = ones(length(f),1) * [c0 cL c0];
rho   = [rho0 rhoS rho0];
[R,T] = slabVibro(f,rho,c,e);
toc

% Comparison (db)
ref = 20*log10(abs([R ; T])); 
sol = 20*log10(abs(Pmes));

norm(ref-sol)/norm(ref)



disp('~~> Michto gypsilab !')



