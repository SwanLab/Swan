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
%|    #    |   FILE       : nrtBmmMaxwellCFIE.m                           |
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
tol  = 1e-3
Ninf = 1e3
beta = 0.8

% Incident direction and field
X0 = [0 0 -1]; 
E  = [0 1  0]; % Polarization (+x for Theta-Theta and +y for Phi-Phi)
H  = cross(X0,E);

% Spherical mesh
mesh = msh('../../meshes/unitSphere_1e3.ply');

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

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k) ;
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k) ;
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k) ;

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
    for j = 1:Nlab   
        % Mass matrix
        Id = integral(sigmaLab{labindex},uLab{labindex},uLab{j});
        
        % Single layer
        T = 1i*k/(4*pi) .* integral( sigmaLab{labindex} , sigmaLab{j} , ...
            uLab{labindex} , Gxy , uLab{j} , tol) ...
            - 1i/(4*pi*k) .* integral( sigmaLab{labindex} ,  sigmaLab{j} , ...
            div(uLab{labindex}) , Gxy , div(uLab{j}) , tol) ;
        Tr = 1i*k/(4*pi) .* regularize( sigmaLab{labindex} , sigmaLab{j} , ...
            uLab{labindex} , '[1/r]' , uLab{j}) ...
            - 1i/(4*pi*k) .* regularize( sigmaLab{labindex} , sigmaLab{j} , ...
            div(uLab{labindex}) , '[1/r]' , div(uLab{j}));
        
        % Double layer
        nxK  = 1/(4*pi) .* integral( sigmaLab{labindex} , sigmaLab{j} , ...
            nx(uLab{labindex}) , gradyGxy , uLab{j} , tol);
        nxKr = 1/(4*pi) .* regularize( sigmaLab{labindex} , sigmaLab{j} , ...
            nx(uLab{labindex}) , 'grady[1/r]' , uLab{j});
        
        % Left hand side
        M{j} = - beta.*(T+Tr)  + (1-beta).*(0.5*Id-(nxK+nxKr));
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
P = bmm(Ilab,Ilab,M);
toc

% Factorization for preconditionning
tic
[L,U] = lu(P);
toc

% Graphical representation
figure
subplot(2,2,1:2)
spy(P)
subplot(2,2,3)
spy(L)
subplot(2,2,4)
spy(U)

% Solve linear system with preconditionner
tic
J = mgcr(LHS,RHS,[],tol,100,L,U);
toc


%%% SURFACIC RADIATION
disp('~~~~~~~~~~~~~ SURFACIC RADIATION ~~~~~~~~~~~~~')

% Mesh Interpolation
Jmsh = feval(u,double(J),mesh);
V    = sqrt(sum(real(cell2mat(Jmsh)).^2,2));

% Graphical representation
figure
% plot(mesh)
hold on
plot(mesh,V)
axis equal
title('|J| surfacic')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar
% caxis([0 3])
view(20,30)
camlight
material dull


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/Ninf .* (1:Ninf)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Tinf = integral(nu,sigma,Ginf,u);
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



disp('~~> Michto gypsilab !')
