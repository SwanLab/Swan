function I = domIntegral3(data)
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
%|    #    |   FILE       : domIntegral3.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 3 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% FUNCTION INTEGRATION --> \int_{domain(y)} f(x,y) dy
if isnumeric(data{1}) && isa(data{2},'dom')
    % Evaluation points
    X  = data{1};
    Nx = size(X,1);
    
    % Domain with quadrature
    Ydom   = data{2};
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    
    % Function evaluation on quadrature
    F = data{3};
    if iscell(F)
        Fxy = {zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)};
        for j = 1:Ny
            Fxy{1}(:,j) = F{1}(X,Y(j,:));
            Fxy{2}(:,j) = F{2}(X,Y(j,:));
            Fxy{3}(:,j) = F{3}(X,Y(j,:));
        end
    else
        Fxy = zeros(Nx,Ny);
        for j = 1:Ny
            Fxy(:,j) = F(X,Y(j,:));
        end
    end
    
    % Integration
    if iscell(Fxy)
        I{1} = Fxy{1} * Wy;
        I{2} = Fxy{2} * Wy;
        I{3} = Fxy{3} * Wy;
    else
        I = Fxy * Wy;
    end

    
%%% FUNCTION INTEGRATION --> \int_{domain(x)} f(x,y) dx
elseif isa(data{1},'dom') && isnumeric(data{2})
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    
    % Evaluation points
    Y  = data{2};
    Ny = size(Y,1);
    
    % Function evaluation on quadrature
    F = data{3};
    if iscell(F)
        Fxy = {zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny)};
        for j = 1:Ny
            Fxy{1}(:,j) = F{1}(X,Y(j,:));
            Fxy{2}(:,j) = F{2}(X,Y(j,:));
            Fxy{3}(:,j) = F{3}(X,Y(j,:));
        end
    else
        Fxy = zeros(Nx,Ny);
        for j = 1:Ny
            Fxy(:,j) = F(X,Y(j,:));
        end
    end
    
    % Integration
    if iscell(Fxy)
        I{1} = Wx' * Fxy{1};
        I{2} = Wx' * Fxy{2};
        I{3} = Wx' * Fxy{3};
    else
        I = Wx' * Fxy;
    end
    
    
%%% FINITE ELEMENT INTEGRATION --> \int_{domain(x)} f(x) psi(x) dx
elseif isa(data{1},'dom') && (isa(data{2},'function_handle') || iscell(data{2})) && isa(data{3},'fem')
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Function evaluation on quadrature
    F = data{2};
    if iscell(F)
        Fx{1} = F{1}(X).';
        Fx{2} = F{2}(X).';
        Fx{3} = F{3}(X).';
    else
        Fx = F(X).';
    end
    
    % Integrated finite element matrix
    v  = data{3};
    Mv = v.uqm(Xdom);
    if iscell(Mv)
        Mv{1} = Wx * Mv{1};
        Mv{2} = Wx * Mv{2};
        Mv{3} = Wx * Mv{3};
    else
        Mv = Wx * Mv;
    end
    
    % Integration
    I = femMultiplyCell(Fx,Mv);
    
    
%%% FINITE ELEMENT INTEGRATION --> \int_{domain(x)} psi(x)' f(x)  dx
elseif isa(data{1},'dom') && isa(data{2},'fem') && (isa(data{3},'function_handle') || iscell(data{3}))
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
        Fx{1} = F{1}(X);
        Fx{2} = F{2}(X);
        Fx{3} = F{3}(X);
    else
        Fx = F(X);
    end
    
    % Integration
    I = femMultiplyCell(Mu,Fx);
    
    
%%% FINITE ELEMENT OPERATOR --> \int_{domain(x)} psi'(x) psi(x) dx
elseif isa(data{1},'dom') && isa(data{2},'fem') && isa(data{3},'fem')
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
    
    % Finite element matrix
    v  = data{3};
    Mv = v.uqm(Xdom);
    
    % Integration
    I = femMultiplyCell(Mu,Mv);
    
    
else
    error('domIntegral3.m : unavailable case')
end
end
