function [Rm1,rRm1,gradRm1,gradrRm1] = domSemiAnalyticInt3D(X0,S,n,tau,nu) 
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : domSemiAnalytic3D.m                           |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2019                                    |
%| ( === ) |   SYNOPSIS   : Semi-analytic integration on triangle for a   |
%|  `---'  |                set of particles                              |
%+========================================================================+

% Vector particles X0 -> Vertices of the triangle
NX0   = size(X0,1);
unNX0 = ones(NX0,1);
un3   = ones(1,3);
xX0S  = unNX0*(S(:,1)') - X0(:,1)*un3;
yX0S  = unNX0*(S(:,2)') - X0(:,2)*un3;
zX0S  = unNX0*(S(:,3)') - X0(:,3)*un3;

% Vector norm
nrmX0S = sqrt(xX0S.^2 + yX0S.^2 + zX0S.^2);

% Height of particles to triangle
h = [xX0S(:,1),yX0S(:,1),zX0S(:,1)]*n';

% Soli angle (cf wikipedia)
ca = zeros(NX0,3);
for i = 1:3
    ip1 = mod(i,3) + 1;
    ca(:,i) = 1./(nrmX0S(:,i).*nrmX0S(:,ip1)) .* ( ...
        xX0S(:,i).*xX0S(:,ip1) + yX0S(:,i).*yX0S(:,ip1) + zX0S(:,i).*zX0S(:,ip1) ) ;
end
sa = sqrt(1-ca.^2);
omega = -pi;
for i = 1:3
    ip1 = mod(i,3) + 1;
    ip2 = mod(ip1,3) + 1;
    omega = omega + acos( (ca(:,i)-ca(:,ip1).*ca(:,ip2)) ./ (sa(:,ip1).*sa(:,ip2)) );
end
omega = omega.*sign(h);
omega(abs(h)<=1e-8) = 0;

% Output initialization
Rm1      = -h.*omega;
rRm1     = zeros(NX0,3);
gradRm1  = -omega*n;
gradrRm1 = zeros(NX0,3,3);

% Edge integration
for ia = 1:3
    % Edge numerotation
    iap1 = mod(ia,3)+1;
    iap2 = mod(iap1,3)+1;
    
    % Projections
    ps   = xX0S(:,iap1)*tau(ia,1) + yX0S(:,iap1)*tau(ia,2) + zX0S(:,iap1)*tau(ia,3);
    psp1 = xX0S(:,iap2)*tau(ia,1) + yX0S(:,iap2)*tau(ia,2) + zX0S(:,iap2)*tau(ia,3);
    ps2  = ps .* (nrmX0S(:,iap2)-nrmX0S(:,iap1));
    
    % Length and height of the edge
    ar = abs(psp1(1)-ps(1));
    xnu = [(xX0S(:,iap1) - ps.*tau(ia,1)), ...
           (yX0S(:,iap1) - ps.*tau(ia,2)), ...
           (zX0S(:,iap1) - ps.*tau(ia,3))] ;
    ah = sqrt(sum(xnu.^2,2));
    
    % Integration of 1/r on the edge - general case h>1e-8
    intaRm1 = asinh(psp1./ah) - asinh(ps./ah);
    
    % Integration of 1/r on the edge - singular case 
    im = (nrmX0S(:,iap1)-abs(ps) < 1e-8);
    intaRm1(im) = 0;
    il = logical((ps>1e-8) .* (psp1>1e-8) .* im);    % particles on the left of the edge
    intaRm1(il) = log(nrmX0S(il,iap2)./nrmX0S(il,iap1));
    ir = logical((ps<-1e-8) .* (psp1<-1e-8) .* im);  % paticles of the right of the edge
    intaRm1(ir) = -log(nrmX0S(ir,iap2)./nrmX0S(ir,iap1)); 
    
    % terme Cn (terme dependant de la position du point d'observation / arete)
    inta = (sqrt(ah.^2 + psp1.^2)-sqrt(ah.^2 + ps.^2));
    inta_rRm1(:,1) = xnu(:,1) .* intaRm1 + inta * tau(ia,1);
    inta_rRm1(:,2) = xnu(:,2) .* intaRm1 + inta * tau(ia,2);
    inta_rRm1(:,3) = xnu(:,3) .* intaRm1 + inta * tau(ia,3);

    % Integration r on the edge
    intaR = 0.5.*(ar.*nrmX0S(:,iap2) + intaRm1.*ah.^2 + ps2);
    
    % Integration of 1/r
    Rm1 = Rm1 + intaRm1 .* (...
        xX0S(:,iap1)*nu(ia,1) + ...
        yX0S(:,iap1)*nu(ia,2) + ...
        zX0S(:,iap1)*nu(ia,3) );
    
    % Integration of r/r       
    rRm1 = rRm1 + intaR*nu(ia,:);
    
    % Integration of grad(1/|r|)
    gradRm1 = gradRm1 + intaRm1*nu(ia,:);    
    
    % Integration of grad(r/|r|)
    gradrRm1(:,:,1) = gradrRm1(:,:,1) + nu(ia,1) * inta_rRm1(:,:);
    gradrRm1(:,:,2) = gradrRm1(:,:,2) + nu(ia,2) * inta_rRm1(:,:);
    gradrRm1(:,:,3) = gradrRm1(:,:,3) + nu(ia,3) * inta_rRm1(:,:);
end

% Normal component of grad(r)
rRm1 = rRm1 + (h.*Rm1)*n;
for i = 1:3
    for j = 1:3
        gradrRm1(:,i,j) = gradrRm1(:,i,j) + n(j)*(Rm1(:)*n(i)) + n(j)*h.*gradRm1(:,i);
    end
end
end
