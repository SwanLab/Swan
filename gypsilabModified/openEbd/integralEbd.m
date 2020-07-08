function [I,loc] = integralEbd(Xdom,Ydom,u,green,k,v,tol)
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
%|    #    |   FILE       : integralEbd.m                                 |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Martin Averseng             |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 8 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% EFFCIENT BESSEL DECOMPOSITION WITH BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
% Domain with quadrature
[X,Wx] = Xdom.qud;
Nx     = size(X,1);
Wx     = spdiags(Wx,0,Nx,Nx);

% Domain with quadrature
[Y,Wy] = Ydom.qud;
Ny     = size(Y,1);
Wy     = spdiags(Wy,0,Ny,Ny);

% Finite element matrix with integration
Mu = u.uqm(Xdom);
if iscell(Mu)
    Mu{1} = Mu{1}' * Wx;
    Mu{2} = Mu{2}' * Wx;
    Mu{3} = Mu{3}' * Wx;
else
    Mu = Mu' * Wx;
end

% Finite element matrix with integration
Mv = v.uqm(Ydom);
if iscell(Mv)
    Mv{1} = Wy * Mv{1};
    Mv{2} = Wy * Mv{2};
    Mv{3} = Wy * Mv{3};
else
    Mv = Wy * Mv;
end

% EBD matrix-vector product
if ~iscell(Mu) && ~iscell(green) && ~iscell(Mv)
    [Gxy,loc] = MVproduct(green,X(:,1:2),Y(:,1:2),tol,k);
    I         = @(V) Mu*Gxy(Mv*V);
    loc       = Mu*loc*Mv;
    
% elseif iscell(Mu) && ~iscell(green) && iscell(Mv)
%     I = 0;
%     for i = 1:3
%         I = I + Mu{i} * ffmProduct(X,Y,Mv{i}*V,green,k,tol);
%     end
%     
% elseif iscell(Mu) && iscell(green) && ~iscell(Mv)
%     I = 0;
%     for i = 1:3
%         I = I + Mu{i} * ffmProduct(X,Y,Mv*V,green{i},k,tol);
%     end
%     
% elseif ~iscell(Mu) && iscell(green) && iscell(Mv)
%     I = 0;
%     for i = 1:3
%         I = I + Mu * ffmProduct(X,Y,Mv{i}*V,green{i},k,tol);
%     end
%     
% elseif iscell(Mu) && iscell(green) && iscell(Mv)
%     I   = 0;
%     ind = [1 2 3 ; 1 3 2 ; 2 3 1 ; 2 1 3 ; 3 1 2 ; 3 2 1];
%     sgn = [+1 -1 +1 -1 +1 -1];
%     for i = 1:6
%         I = I + sgn(i) .* (Mu{ind(i,1)} * ...
%             ffmProduct(X,Y,Mv{ind(i,3)}*V,green{ind(i,2)},k,tol));
%     end
    
else
    error('integralEbd.m : unavailable case')
end
end


function MV = test()


end


