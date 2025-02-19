function I = domIntegral4(data)
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
%|    #    |   FILE       : domIntegral4.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 4 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% FINITE ELEMENT OPERATOR --> \int_{mesh(x)} psi(x)' f(x) psi(x) dx 
if isa(data{1},'dom') && isa(data{2},'fem')
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Integrated finite element matrix
    u  = data{2};
    Mu = u.uqm(Xdom);
    if iscell(Mu)
        Mu{1} = Mu{1}' * Wx;
        Mu{2} = Mu{2}' * Wx;
        Mu{3} = Mu{3}' * Wx;
    else
        Mu = Mu' * Wx;
    end
    
    % Function evaluation on quadrature
    F = data{3};
    if iscell(F)
        Fx{1} = spdiags(F{1}(X),0,Nx,Nx);
        Fx{2} = spdiags(F{2}(X),0,Nx,Nx);
        Fx{3} = spdiags(F{3}(X),0,Nx,Nx);
    else
        Fx = spdiags(F(X),0,Nx,Nx);
    end
    
    % Finite element matrix
    v  = data{4};
    Mv = v.uqm(Xdom);
    
    % Integration
    I = femMultiplyCell(Mu,Fx,Mv);

    
%%% BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(y)} f(x,y) psi(y) dy
elseif isnumeric(data{1}) && isa(data{2},'dom')
    % Evaluation points
    X  = data{1};
    Nx = size(X,1);
    
    % Domain with quadrature
    Ydom   = data{2};
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);

    % Function evaluation on quadrature
    F = data{3};
    if iscell(F)
        Fxy = {zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)};
        for i = 1:3
            for j = 1:Ny
                Fxy{i}(:,j) = F{i}(X,Y(j,:));
            end
        end
    else
        Fxy = zeros(Nx,Ny);
        for j = 1:Ny
            Fxy(:,j) = F(X,Y(j,:));
        end
    end    

    % Integrated finite element matrix
    v  = data{4};
    Mv = v.uqm(Ydom);
    if iscell(Mv)
        Mv{1} = Wy * Mv{1};
        Mv{2} = Wy * Mv{2};
        Mv{3} = Wy * Mv{3};
    else
        Mv = Wy * Mv;
    end
        
    % Integration
    I = femMultiplyCell(Fxy,Mv);

    
%%% BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x,y) dx
elseif isa(data{1},'dom') && isnumeric(data{2})
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Evaluation points
    Y  = data{2};
    Ny = size(Y,1);
    
    % Integrated finite element matrix
    u  = data{3};
    Mu = u.uqm(Xdom);
    if iscell(Mu)
        Mu{1} = Mu{1}' * Wx;
        Mu{2} = Mu{2}' * Wx;
        Mu{3} = Mu{3}' * Wx;
    else
        Mu = Mu' * Wx;
    end
    
    % Function evaluation on quadrature
    F = data{4};
    if iscell(F)
        Fxy = {zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)};
        for i = 1:3
            for j = 1:Ny
                Fxy{i}(:,j) = F{i}(X,Y(j,:));
            end
        end
    else
        Fxy = zeros(Nx,Ny);
        for j = 1:Ny
            Fxy(:,j) = F(X,Y(j,:));
        end
    end
    
    % Integration
    I = femMultiplyCell(Mu,Fxy);
    
    
else
    error('domIntegral4.m : unavailable case')
end
end
