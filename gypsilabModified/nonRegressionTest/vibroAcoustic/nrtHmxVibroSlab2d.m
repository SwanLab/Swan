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
%|    #    |   FILE       : nrtHmxVibroSlab2d.m                           |
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

% Accuracy
tol = 1e-3

% Width of the slab
e = 0.5

% Frequency
f = 900:500:4000

% Incident direction (from bottom)
X0  = [0 1 0];

% Exterior domain (water)
rho0 = 1000;         % density      (kg.m3)
c0   = 1500;         % celerity     (m.s-1)
k0   = 2*pi/c0.*f;   % wave-number  (m-1)
lam0 = c0./f;        % wave-length  (m)
    
% Interior domain (different from water)
rhoS = 2*rho0;                                 % density                             (kg.m3)
cL   = 2*c0;                                   % celerity of longitudinal waves      (m.s-1)
kL   = 2*pi.*f./cL;                            % wave-number of longitudinal waves   (m-1)
lamL = real(cL)./f;                            % wavelength of longitudinal waves    (m)
cT   = 0;                                      % celerity of transverse waves        (m.s-1)
kT   = 2*pi.*f./cT;                            % wave-number of transverse waves     (m-1)
lamT = real(cT)./f;                            % wavelength of transverse waves      (m)

% Solution (pressure)
sol = zeros(2,length(f)); 

% Loop for each frequency
for i = 1:length(f)
    % Minimum wavelength 
    tmp  = [lam0(i),lamL(i),lamT(i)];
    lmin = min(tmp(tmp>0));
    
    % Slab mesh
    L    = 60 * lmin;            % 60 wavelength to simulate infinite slab
    nx   = ceil(L/lmin * 6)+1;   % 6 node per wavelength for L
    ny   = ceil(e/lmin * 12)+1;  % 12 node per wavelength for e
    N    = nx * ny;              % Total number of nodes
    mesh = mshSquare(N,[L e])

    % Radiative mesh (fixed number of nodes)
    radiat = mshSquare(1e3,[L L]);
    
    % Boundary
    bound = swap(mesh.bnd)
    
    % Measurement points for trans and refl coeff (1 wavelenth from bound)
    Xmes = [0 -e/2-lmin 0 ; 0 e/2+lmin 0];
    
    % Cut-off function (50% full, 10% decrease) 
    cutoff = vibsCutoff(1,L/5,L/10);
    
    % Green kernel function
    Gxy         = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k0(i));
    gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k0(i));
    gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k0(i));
    gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k0(i));

    % Plane wave function
    PW         = @(X) exp(1i*k0(i)*X*X0') .* cutoff(X);
    gradxPW{1} = @(X) 1i*k0(i)*X0(1) .* PW(X);
    gradxPW{2} = @(X) 1i*k0(i)*X0(2) .* PW(X);
    gradxPW{3} = @(X) 1i*k0(i)*X0(3) .* PW(X);
    
    % Coupling coeff for Brackage-Werner simulation  
    beta = 1i*k0(i);
    
    % Quadrature and finite elements (volumn)
    omega = dom(mesh,3);
    U     = fem(mesh,'P1');

    % Quadrature and finite elements (boundary)
    sigma = dom(bound,3);
    u     = fem(bound,'P1');
    
    % Left-hand side
    tic
    [A,B,C,D] = vibsHmxBlockOperator(omega,U,sigma,u,cL,cT,rhoS,c0,rho0,f(i),tol);
    toc
    
    % Add dirichlet condition to x unknows (penalization)
    A(sub2ind(size(A),1:length(U),1:length(U))) = 1e15;
    
    % Right-hand side
    V    = cell(3,1);
    V{1} = - integral(sigma,ntimes(U,1),PW);
    V{2} = - integral(sigma,ntimes(U,2),PW);
    V{3} = integral(sigma,ntimes(u),gradxPW);
    
    % Resolution with Schur complement
    Fa     = decomposition(A);
    LHS    = @(V) D*V - C*(Fa \ (B.Ml*(B.Mr*V)) );      
    RHS    = V{end} - C*(Fa \ cell2mat(V(1:end-1)) );
    mu     = mgcr(LHS,RHS,[],tol,100); 
    lambda = beta*mu;
    
    % Measure of refexive and transmitted coeff
    tic
    Pmes = 1i/4 .* integral(Xmes,sigma,Gxy,u)*lambda - ...
        1i/4 .* integral(Xmes,sigma,gradyGxy,ntimes(u))*mu;
    Pmes(2) = Pmes(2) + PW(Xmes(2,:));
    toc
    
    % Save solution
    sol(:,i) = Pmes;
end

% Analytical solution
tic
c     = ones(length(f),1) * [c0 cL c0];
rho   = [rho0 rhoS rho0];
[R,T] = slabVibro(f,rho,c,e);
toc

% Comparison (db)
ref = 20*log10(abs([R ; T])); 
sol = 20*log10(abs(sol));

norm(ref-sol)/norm(ref)

% Graphical representation
figure(100)
subplot(1,2,1)
plot(f,ref(1,:),'r',f,sol(1,:),'b+')
grid on
title('Reflexion coeffiscient')
legend({'Analytical','Numerical'})
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

subplot(1,2,2)
plot(f,ref(2,:),'r',f,sol(2,:),'b+')
grid on
title('Transmission coeffiscient')
legend({'Analytical','Numerical'})
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')




disp('~~> Michto gypsilab !')


