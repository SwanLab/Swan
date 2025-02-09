%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : nrtOprMaxwellCFIE.m                           |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Solve PEC scatering problem with              |
%|  `---'  |                CFIE formulation                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
run('../../addpathGypsilab.m')

% Parameters
N    = 1e3;
tol  = 1e-3
beta = 0.5

% Incident direction and field
X0 = [0 0 -1]; 
E  = [0 1  0]; % Polarization (+x for Theta-Theta and +y for Phi-Phi)
H  = cross(X0,E);

% Spherical mesh
mesh = mshSphere(N,1);

% Graphical representation
figure
plot(mesh)
axis equal

% Frequency adjusted to maximum esge size
stp = mesh.stp;
k   = 1/stp(2);
c   = 299792458;
f   = (k*c)/(2*pi);
disp(['Frequency : ',num2str(f/1e6),' MHz']);

% Domain
sigma = dom(mesh,3);

% Finite elements
u = fem(mesh,'RWG');

% Incident Plane wave (electromagnetic field)
PWE{1} = @(X) exp(1i*k*X*X0') * E(1);
PWE{2} = @(X) exp(1i*k*X*X0') * E(2);
PWE{3} = @(X) exp(1i*k*X*X0') * E(3);

PWH{1} = @(X) exp(1i*k*X*X0') * H(1);
PWH{2} = @(X) exp(1i*k*X*X0') * H(2);
PWH{3} = @(X) exp(1i*k*X*X0') * H(3);

% Incident wave representation
plot(mesh,real(PWE{2}(mesh.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
view(0,10)


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATORS ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy      = '[exp(ikr)/r]';
gradxGxy = {'gradx[exp(ikr)/r]1','gradx[exp(ikr)/r]2','gradx[exp(ikr)/r]3'};
gradyGxy = {'grady[exp(ikr)/r]1','grady[exp(ikr)/r]2','grady[exp(ikr)/r]3'};

% Number of pool
Nlab = length(Composite)

% Domain decomposition
tic
[Ilab,sigmaLab,uLab] = femSubdivide(sigma,u,Nlab,10);
toc

% Parallel loop 
tic
spmd
    M = cell(1,numlabs);
    P = cell(1,numlabs);
    for j = 1:Nlab   
        % Mass matrix
        Id = integral(sigmaLab{labindex},uLab{labindex},uLab{j});
        
        % EFIE
        T  = oprIntegral('T',k,sigmaLab{labindex},uLab{labindex},sigmaLab{j},uLab{j},tol);
        Tr = 1i*k/(4*pi).*regularize(sigmaLab{labindex},sigmaLab{j},uLab{labindex},'[1/r]',uLab{j}) ...
            - 1i/(4*pi*k).*regularize(sigmaLab{labindex},sigmaLab{j},div(uLab{labindex}),'[1/r]',div(uLab{j}));

        % MFIE
        nxK  = oprIntegral('nxK',k,sigmaLab{labindex},uLab{labindex},sigmaLab{j},uLab{j},tol);
        nxKr = 1/(4*pi).*regularize(sigmaLab{labindex},sigmaLab{j},nx(uLab{labindex}),'grady[1/r]',uLab{j});
        
        % Left hand side
        M{j} = - beta.*(T+Tr)  + (1-beta).*(0.5*Id-(nxK+nxKr));
        P{j} = sparse(M{j},Id);
    end
end
toc

% Define LHS
LHS = @(V) spmdProduct(Ilab,M,V);

% Right hand side
RHS = beta*integral(sigma,u,PWE) - (1-beta)*integral(sigma,nx(u),PWH);


%%% SOLVE LINEAR SYSTEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR SYSTEM ~~~~~~~~~~~~~')

% Block matrix
tic
P = sparse(bmm(Ilab,Ilab,P));
toc

% Factorization for preconditionning
tic
[L,U] = ilu(P);
toc

% Solve linear system with preconditionner
tic
J = mgcr(LHS,RHS,[],tol,100,L,U);
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
Ninf  = 1e3;
theta = 2*pi/1e3 .* (1:Ninf)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
Ginf = '[exp(-ikxy)]';

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Tinf = integral(nu,sigma,Ginf,k,u,tol);
sol  = 1i*k/(4*pi)*cross(nu, cross([Tinf{1}*J, Tinf{2}*J, Tinf{3}*J], nu));

% Radiation infinie de reference, convention e^(+ikr)/r
nMax = 100; refInf = zeros(Ninf,1);
if E(1) == 1
    for jj = 1:Ninf
        refInf(jj,:) = sphereMaxwell(1, -f, theta(jj), 0.0, nMax);
    end
else
    for jj = 1:Ninf
        [~,refInf(jj)] = sphereMaxwell(1, -f, theta(jj), pi/2, nMax);
    end
end
refInf = refInf ./ sqrt(4*pi);

% Radiations infinies en X
if E(1) == 1
    sol = sin(theta)'.*sol(:,3) - cos(theta)'.*sol(:,1);
else
    sol = sol(:,2);
end

% Erreur
eL2   = norm(refInf-sol,2)/norm(refInf,2)
eLINF = norm(refInf-sol,'inf')/norm(refInf,'inf')
             
% Representation graphique
figure
subplot(1,2,1)
plot(theta,20*log10(abs(sol)),'b',theta,20*log10(abs(refInf)),'r--')

subplot(1,2,2)
plot(theta,real(sol),'--b', theta,imag(sol),'--r', theta, real(refInf),':b', theta,imag(refInf),':r');
drawnow


%%% SURFACIC RADIATION
disp('~~~~~~~~~~~~~ SURFACIC RADIATION ~~~~~~~~~~~~~')

% Mesh Interpolation
Jmsh = feval(u,J,mesh);
V    = sqrt(sum(real(cell2mat(Jmsh)).^2,2));

% Graphical representation
figure
plot(mesh)
hold on
plot(mesh,V)
axis equal
title('|J| surfacic')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar



disp('~~> Michto gypsilab !')


