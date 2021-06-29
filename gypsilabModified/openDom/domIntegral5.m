function I = domIntegral5(data)
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
%|    #    |   FILE       : domIntegral5.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 5 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
if isa(data{1},'dom') && isa(data{2},'dom')
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Domain with quadrature
    Ydom   = data{2};
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);
    
    % Finite element matrix with integration
    u  = data{3};
    Mu = u.uqm(Xdom);
    if iscell(Mu)
        Mu{1} = Mu{1}' * Wx;
        Mu{2} = Mu{2}' * Wx;
        Mu{3} = Mu{3}' * Wx;
    else
        Mu = Mu' * Wx;
    end  
    
    % Function applied to quadrature
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
    
    % Finite element matrix with integration
    v  = data{5};
    Mv = v.uqm(Ydom);
    if iscell(Mv)
        Mv{1} = Wy * Mv{1};
        Mv{2} = Wy * Mv{2};
        Mv{3} = Wy * Mv{3};
    else
        Mv = Wy * Mv;
    end
    
    % Integration
    I = femMultiplyCell(Mu,Fxy,Mv);
    
    
%%% H-MATRIX BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(y)} f(x,y) psi(y) dy
elseif isnumeric(data{1}) && isa(data{2},'dom')
    % Evaluation points 
    X  = data{1};
    Nx = size(X,1);
    Mx = speye(Nx,Nx);
    
    % Domain with quadrature
    Ydom   = data{2};
    [Y,Wy] = Ydom.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);
    
    % Green kernel
    green = data{3};
    
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
    
    % Accuracy
    tol = data{5};
        
    % H-Matrix integration
    if iscell(Mv) && ~iscell(green)
        I{1} = hmx(X,v.unk,Mx,X,green,Y,Mv{1},tol);
        I{2} = hmx(X,v.unk,Mx,X,green,Y,Mv{2},tol);
        I{3} = hmx(X,v.unk,Mx,X,green,Y,Mv{3},tol);
        
    elseif ~iscell(Mv) && iscell(green)
        I{1} = hmx(X,v.unk,Mx,X,green{1},Y,Mv,tol);
        I{2} = hmx(X,v.unk,Mx,X,green{2},Y,Mv,tol);
        I{3} = hmx(X,v.unk,Mx,X,green{3},Y,Mv,tol);
        
    else
        I = hmx(X,v.unk,Mx,X,green,Y,Mv,tol);
    end

    
%%% H-MATRIX BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x,y) dx
elseif isa(data{1},'dom') && isnumeric(data{2})
    % Domain with quadrature
    Xdom   = data{1};
    [X,Wx] = Xdom.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Evaluation points
    Y  = data{2};
    Ny = size(Y,1);
    My = speye(Ny,Ny);
    
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
    
    % Green kernel
    green = data{4};
    
    % Accuracy
    tol = data{5};
    
    % H-Matrix integration
    if iscell(Mu) && ~iscell(green)
        I{1} = hmx(u.unk,Y,Mu{1},X,green,Y,My,tol);
        I{2} = hmx(u.unk,Y,Mu{2},X,green,Y,My,tol);
        I{3} = hmx(u.unk,Y,Mu{3},X,green,Y,My,tol);
        
    elseif ~iscell(Mu) && iscell(green)
        I{1} = hmx(u.unk,Y,Mu,X,green{1},Y,My,tol);
        I{2} = hmx(u.unk,Y,Mu,X,green{2},Y,My,tol);
        I{3} = hmx(u.unk,Y,Mu,X,green{3},Y,My,tol);
        
    else
        I = hmx(u.unk,Y,Mu,X,green,Y,My,tol);
    end

else
    error('domIntegral5.m : unavailable case')
end
end
