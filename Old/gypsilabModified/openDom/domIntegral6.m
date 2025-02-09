function I = domIntegral6(data)
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
%|    #    |   FILE       : domIntegral6.m                                |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 6 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-MATRIX BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
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
    
    % Green kernel
    green = data{4};
    
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
    
    % Accuracy
    tol = data{6};
    
    % H-Matrix Integration
    if iscell(Mu) && ~iscell(green) && ~iscell(Mv)
        I{1} = hmx(u.unk,v.unk,Mu{1},X,green,Y,Mv,tol);
        I{2} = hmx(u.unk,v.unk,Mu{2},X,green,Y,Mv,tol);
        I{3} = hmx(u.unk,v.unk,Mu{3},X,green,Y,Mv,tol);
        
    elseif ~iscell(Mu) && iscell(green) && ~iscell(Mv)
        I{1} = hmx(u.unk,v.unk,Mu,X,green{1},Y,Mv,tol);
        I{2} = hmx(u.unk,v.unk,Mu,X,green{2},Y,Mv,tol);
        I{3} = hmx(u.unk,v.unk,Mu,X,green{3},Y,Mv,tol);
        
    elseif ~iscell(Mu) && ~iscell(green) && iscell(Mv)
        I{1} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{1},tol);
        I{2} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{2},tol);
        I{3} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{3},tol);
        
    else
        I = hmx(u.unk,v.unk,Mu,X,green,Y,Mv,tol);
    end


%%% FFM BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(y)} f(x,y) psi(y) dy
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
    
    % Wave number
    k = data{4};
    
    % Integrated finite element matrix
    v  = data{5};
    Mv = v.uqm(Ydom);
    if iscell(Mv)
        Mv{1} = Wy * Mv{1};
        Mv{2} = Wy * Mv{2};
        Mv{3} = Wy * Mv{3};
    else
        Mv = Wy * Mv;
    end
    
    % Accuracy
    tol = data{6};
        
    % Fast & Furious integration
    if iscell(Mv) && iscell(green)
        I = ffm(Mx,X,green{1},k,Y,Mv{1},tol) + ...
            ffm(Mx,X,green{2},k,Y,Mv{2},tol) + ...
            ffm(Mx,X,green{3},k,Y,Mv{3},tol);
        
    elseif iscell(Mv) && ~iscell(green)
        I{1} = ffm(Mx,X,green,k,Y,Mv{1},tol);
        I{2} = ffm(Mx,X,green,k,Y,Mv{2},tol);
        I{3} = ffm(Mx,X,green,k,Y,Mv{3},tol);
        
    elseif ~iscell(Mv) && iscell(green)
        I{1} = ffm(Mx,X,green{1},k,Y,Mv,tol);
        I{2} = ffm(Mx,X,green{2},k,Y,Mv,tol);
        I{3} = ffm(Mx,X,green{3},k,Y,Mv,tol);
        
    else
        I = ffm(Mx,X,green,k,Y,Mv,tol);
    end

    
%%% FFM BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x,y) dx
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
    
    % Wave number
    k = data{5};

    % Accuracy
    tol = data{6};
    
    % Fast & Furious integration
    if iscell(Mu) && iscell(green)
        I = ffm(Mu{1},X,green{1},k,Y,My,tol) + ...
            ffm(Mu{2},X,green{2},k,Y,My,tol) + ...
            ffm(Mu{3},X,green{3},k,Y,My,tol) ;
        
    elseif iscell(Mu) && ~iscell(green)
        I{1} = ffm(Mu{1},X,green,k,Y,My,tol);
        I{2} = ffm(Mu{2},X,green,k,Y,My,tol);
        I{3} = ffm(Mu{3},X,green,k,Y,My,tol);
        
    elseif ~iscell(Mu) && iscell(green)
        I{1} = ffm(Mu,X,green{1},k,Y,My,tol);
        I{2} = ffm(Mu,X,green{2},k,Y,My,tol);
        I{3} = ffm(Mu,X,green{3},k,Y,My,tol);
        
    else
        I = ffm(Mu,X,green,k,Y,My,tol);
    end

    
else
    error('domIntegral6.m : unavailable case')
end
end
