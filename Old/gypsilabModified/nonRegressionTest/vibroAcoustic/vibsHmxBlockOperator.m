function [A,B,C,D] = vibsHmxBlockOperator(omega,U,sigma,u,cL,cT,rhoS,c0,rho0,f,tol)
%+========================================================================+
%|                                                                        |
%|                 OPENVIBS - LIBRARY FOR VIBRO-ACOUSTIC                  |
%|           openVibs is part of the GYPSILAB toolbox for Matlab          |
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
%|    #    |   FILE       : vibsHmxBlockOperator.m                        |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Constants
mu     = rhoS*cT^2;
lambda = rhoS*(cL^2 - 2*cT^2);
w      = 2*pi*f;
k      = w/c0;  

% Dimension 
n = size(omega.msh.elt,2)-1;

% Green kernel function
if (n == 2)
    Gxy         = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);
    gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
    gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
    gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);
    G0          = '[log(r)]';
    gradyG0     = 'grady[log(r)]';    
    cteGxy      = 1i/4;
    cteG0       = -1/(2*pi);
    
elseif (n == 3)
    Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
    gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
    gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
    gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);
    G0          = '[1/r]';
    gradyG0     = 'grady[1/r]';    
    cteGxy      =  1/(4*pi);
    cteG0       =  1/(4*pi);
    
else
    error('vibsNeumannBW.m : unavailable case.')
end
    
% Coupling coeff for Brackage-Werner simulation
beta = 1i*k;

%%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTO-ELASTO (A11) %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
A = cell(n,n);

% Static part
GG = integral(omega,grad(U),grad(U));
for i = 1:n
    for j = 1:n        
        % Operator div(U):div(U)
        DD = integral(omega,grad(U,i),grad(U,j));
        
        % Operator e(U):e(U)
        EE = integral(omega,grad(U,j),grad(U,i));
        if (i==j)
            EE = EE + GG;
        end
        
        % Summation
        A{i,j} = lambda.*DD + mu.*EE;
    end
end

% Dynamic part
if (f ~= 0)
    Id = integral(omega,U,U);
    for i = 1:n
       A{i,i} = A{i,i} - (rhoS*w^2) .* Id; 
    end
end

% Final form (sparse)
A = cell2mat(A);

%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTO-ACOUSTIC (A12) %%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian quadrature
[Xqud,Wx] = sigma.qud;
Wx        = spdiags(Wx,0,length(Wx),length(Wx));

% Collocation mass operator
Id = u.uqm(sigma);

% Collocation boundary operator
S = cteGxy .* integral(Xqud,sigma,Gxy,u,tol) + ...
    cteG0  .* regularize(Xqud,sigma,G0,u);

% Collocation boundary operator
D = cteGxy .* integral(Xqud,sigma,gradyGxy,ntimes(u),tol) + ...
    cteG0  .* regularize(Xqud,sigma,gradyG0,ntimes(u));

% Final operator Brackage-Werner : [1i*k*beta*S - (Id/2 + D)]
B.Mr = beta.*S - (0.5*Id + D);

% Normal trace of the volumn element matrix
nU   = ntimes(U);
nPHI = nU.uqm(sigma);

% Coupling to FEM
B.Ml = cell(n,1);
for i = 1:n
    B.Ml{i} = (nPHI{i}' * Wx);
end
B.Ml = cell2mat(B.Ml);

%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTO-ACOUSTIC (A21) %%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
C = cell(1,n);

% Coupling FEM 
for i = 1:n
    C{i} = (rho0*w^2) .* integral(sigma,u,ntimes(U,i));
end

% Final format (sparse)
C = cell2mat(C);

%%%%%%%%%%%%%%%%%%%%%%%%%% ACOUSTIC-ACOUSTIC (A22) %%%%%%%%%%%%%%%%%%%%%%%%

% Finite element mass matrix
Id = integral(sigma,u,u);

% Finite element boundary operator
H  = cteGxy .* (k^2 * integral(sigma,sigma,ntimes(u),Gxy,ntimes(u),tol) ...
    - integral(sigma,sigma,nxgrad(u),Gxy,nxgrad(u),tol));
Hr = cteG0  .* (k^2 * regularize(sigma,sigma,ntimes(u),G0,ntimes(u)) ...
    - regularize(sigma,sigma,nxgrad(u),G0,nxgrad(u)));

% Finite element boundary operator
D  = cteGxy .* integral(sigma,sigma,u,gradyGxy,ntimes(u),tol);
Dr = cteG0  .* regularize(sigma,sigma,u,gradyG0,ntimes(u));

% Final operator Brackage-Werner : - [1i*k*beta*(-Id/2 + Dt) - H]
D = - (beta.*(-0.5*Id + (D+Dr).') - (H+Hr));
end
