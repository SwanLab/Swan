%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Yosra Boukari & Houssem Haddar (c) 2017-2018.|
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtIpbHelmholtz0.m                            |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Yosra Boukari               |
%|  ( # )  |                & Houssem Haddar                              |
%|  / 0 \  |   CREATION   : 14.03.2017                                    |
%| ( === ) |   LAST MODIF : 14.03.2018                                    |
%|  `---'  |   SYNOPSIS   : Completion data for spherical helmholtz       |
%|         |                scatering                                     |
%|         |   ref1 : A Convergent Data Completion Algorithm Using Surface|
%|         |   Integral Equation, Inverse Problems, IOP Publishing ...    |
%|         |   ref2: Poly cours Terasse p. 194, Calderon operators        |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Wave number and frequency
k = 1
f = (k*340)/(2*pi)

% Data noise
noise = 0%1e-2

% Incident plane wave with normlized direction
X0         = [0 0 -1];
X0         = X0./norm(X0);
PW         = @(X) exp(1i*k*X*X0');
gradxPW{1} = @(X) 1i*k*X0(1) .* PW(X);
gradxPW{2} = @(X) 1i*k*X0(2) .* PW(X);
gradxPW{3} = @(X) 1i*k*X0(3) .* PW(X);

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% grady Green kernel function --> G(x,y) = grady[exp(ik|x-y|)/|x-y|]
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Green kernel function --> G(x,y) = gradx[exp(ik|x-y|)/|x-y|]
dxGxy{1} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]1',k);
dxGxy{2} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]2',k);
dxGxy{3} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]3',k);

% Finite element
gss = 3;
typ = 'P1';
tol = 1e-3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATERING PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fictive mesh for diffraction problem (sphere)
mesh = mshSphere(4e2,1);

% % Fictive mesh for diffraction problem (cube)
% mesh = mshCube(4e2,[1 1 1]);
% mesh = mesh.bnd;

% Verify wave number
stp = mesh.stp;
if k > 1/stp(2)
    error('completionHelmholtz.m : wave number k is too high for mesh')
end

% Quadrature
sigma = dom(mesh,gss);

% Finite element
u = fem(mesh,typ);

% Boundary element operator on fictive domain (Hypersingular)
H = 1/(4*pi) .* (k^2 * integral(sigma,sigma,ntimes(u),Gxy,ntimes(u)) ...
    - integral(sigma,sigma,nxgrad(u),Gxy,nxgrad(u)));
H = H + 1/(4*pi) .* (k^2 * regularize(sigma,sigma,ntimes(u),'[1/r]',ntimes(u)) ...
    - regularize(sigma,sigma,nxgrad(u),'[1/r]',nxgrad(u)));

% Solve neumann problem on fictive boundary : - [H] mu = - dnP0
RHS = - integral(sigma,ntimes(u),gradxPW);
mu  = (-H) \ RHS;

% Radiation on boundary :  p(x) =  - [D + I/2] mu
Id  = integral(sigma,u,u);
D   = 1/(4*pi) .* integral(sigma,sigma,u,dyGxy,ntimes(u));
D   = D + 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(u));
sol = - Id \ ( (D + 0.5*Id) * mu);

% Comparare to analytic solution
ref = sphereHelmholtz('dom','neu',1,k,1.0001*mesh.vtx);
errScatMesh = norm(ref - sol)/norm(ref)

% Graphical representation
figure
plot(mesh)
hold on
plot(mesh,abs(sol+PW(u.dof)))
xlabel('X');   ylabel('Y');   zlabel('Z');
title('Total field')
axis equal
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATERING RADIATION %%%%%%%%%%%%%%%%%%%%%%%%%%
% Interior boundary
meshInt = mshSphere(5e2,1.25);
meshInt = swap(meshInt);

% Exterior boundary
meshExt = mshSphere(1e3,1.5);

% Fusion des maillages
meshTot = union(meshExt,meshInt);

% Verify wave number
stp = mesh.stp;
if k > 1/stp(2)
    error('completionHelmholtz.m : wave number k is too high for mesh')
end

% Quadrature
sigmaTot = dom(meshTot,gss);

% Finite element
uTot = fem(meshTot,typ);

% Mass matrix
IdTot = integral(sigmaTot,uTot,uTot);

% Solution of the scatering problem : p(x) = - [D] mu
p = - 1/(4*pi) .* integral(sigmaTot,sigma,uTot,dyGxy,ntimes(u)) * mu;

% Extract Galerkin
p = IdTot \ p;

% Compare to analytical solution
ref = sphereHelmholtz('dom','neu',1,k,uTot.dof);
errScat = norm(ref-p)/norm(ref)

% Scatered speed : dnp(x) = - [H] mu
dnp = - 1/(4*pi) .* (k^2 * integral(sigmaTot,sigma,ntimes(uTot),Gxy,ntimes(u)) ...
    - integral(sigmaTot,sigma,nxgrad(uTot),Gxy,nxgrad(u))) * mu;

% Extract Galerkin
dnp = IdTot \ dnp;

% % Graphical representation
% figure
% plot(mesh)
% hold on
% plot(mesh,abs(p+PW(uTot.dof)))
% plotNrm(mesh,'r')
% xlabel('X');   ylabel('Y');   zlabel('Z');
% title('Interior radiation')
% axis equal
% alpha(0.5)
% colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPERATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single layer --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy
S = 1/(4*pi) .* integral(sigmaTot,sigmaTot,uTot,Gxy,uTot);
S = S + 1/(4*pi) .* regularize(sigmaTot,sigmaTot,uTot,'[1/r]',uTot);

% Double layer --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy
D = 1/(4*pi) .* integral(sigmaTot,sigmaTot,uTot,dyGxy,ntimes(uTot));
D = D + 1/(4*pi) .* regularize(sigmaTot,sigmaTot,uTot,'grady[1/r]',ntimes(uTot));

% Double layer --> \int_Sx \int_Sy psi(x)' dnx G(x,y) psi(y) dx dy
Dt = D.';

% Hypersingular --> k^2 * \int_Sx \int_Sy n.psi(x) G(x,y) n.psi(y) dx dy
%                   - \int_Sx \int_Sy nxgrad(psi(x)) G(x,y) nxgrad(psi(y)) dx dy
N = 1/(4*pi) .* (k^2 * integral(sigmaTot,sigmaTot,ntimes(uTot),Gxy,ntimes(uTot)) ...
    - integral(sigmaTot,sigmaTot,nxgrad(uTot),Gxy,nxgrad(uTot)));
N = N + 1/(4*pi) .* (k^2 * regularize(sigmaTot,sigmaTot,ntimes(uTot),'[1/r]',ntimes(uTot)) ...
    - regularize(sigmaTot,sigmaTot,nxgrad(uTot),'[1/r]',nxgrad(uTot)));

% Calderon operators
Z = sparse(length(uTot),length(uTot));
I = [IdTot Z ; Z IdTot];
H = [-D S ; -N Dt];

% Calderon operator for interior problem (projector)
A = H + 0.5*I;
errCalderon = norm(A*(I\A)-A,'fro')/norm(A,'fro')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Indices des dof
N    = size(meshTot.vtx,1);
Iext = 1:size(meshExt.vtx,1);
Iint = Iext(end)+1:N;

% Excitation
f = p(Iext);
g = dnp(Iext);

% Solutions
pSol   = p(Iint);
dnpSol = dnp(Iint);

% Indice operateur
Iext = [Iext , N + Iext];
Iint = [Iint , N + Iint];

% Extraction des sous matrices
A11 = A(Iext,Iext);
A12 = A(Iext,Iint);
A21 = A(Iint,Iext);
A22 = A(Iint,Iint);

% Matrices identite
I11 = I(Iext,Iext);
I22 = I(Iint,Iint);

% Linear system
LHS = [A12 ; A22 - I22];

% Right hand side
RHS = [ I11 - A11 ; -A21] * [f;g];

% Solve linear system with gaussian noise
X = LHS \ (RHS .* (1 + noise*(1 + randn(size(RHS,1),1))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TICHONOV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data  = F1;
% datan = B * F1n;
% [U,S,V]=svd(A);
% global s2 delta rhs2
% a1 = length(fe{2});
% rhs=(U')*datan;
% rhsr =rhs(1:a1+a1);
% s2 = diag(S).^2;
% rhs2 = abs(rhsr).^2;
% delta=norm(datan-B*(data))/norm(datan);  
% alpha0=delta*min(s2)/(1-delta);
% alpha1=delta*max(s2)/(1-delta);
% 
% fmor = @(x)sum(((x ./(x+s2)).^2 - delta^2).*rhs2);
% alphaM = fzero(fmor,[alpha0 alpha1]); 
% 
% X2 = V*((diag(S)./(alphaM+s2)).*rhsr);
% errTichonov = norm(X-X2)./norm(X2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare external pression
sol = X(1:end/2);
ref = pSol;
errSolution = norm(ref-sol)/norm(ref)

% Graphical representation
figure(13)
plot(meshInt,abs(sol+PW(meshInt.vtx)))
xlabel('X');   ylabel('Y');   zlabel('Z');
title('Reconstructed total field')
axis equal
colorbar

% Graphical representation
figure(14)
plot(meshInt,abs(sol-ref)./abs(ref))
xlabel('X');   ylabel('Y');   zlabel('Z');
title('Relative error on reconstructed field')
axis equal
colorbar

% Compare normal derivative
sol = X(end/2+1:end);
ref = dnpSol;
errDerivative = norm(ref-sol)/norm(ref)


disp('~~> Michto gypsilab !')



