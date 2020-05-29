function I = domIntegral7(data)
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
%|    #    |   FILE       : domIntegral7.m                                |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 8 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% FAST AND FURIOUS METHOD WITH BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
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

% Wave number
k = data{5};

% Finite element matrix with integration
v  = data{6};
Mv = v.uqm(Ydom);
if iscell(Mv)
    Mv{1} = Wy * Mv{1};
    Mv{2} = Wy * Mv{2};
    Mv{3} = Wy * Mv{3};
else
    Mv = Wy * Mv;
end

% Accuracy
tol = data{7};

% FFM matrix-vector product
if ~iscell(Mu) && ~iscell(green) && ~iscell(Mv)
    I = ffm(Mu,X,green,k,Y,Mv,tol);
    
elseif iscell(Mu) && ~iscell(green) && iscell(Mv)
    I = ffm(Mu{1},X,green,k,Y,Mv{1},tol) + ...
        ffm(Mu{2},X,green,k,Y,Mv{2},tol) + ...
        ffm(Mu{3},X,green,k,Y,Mv{3},tol) ;
    
elseif iscell(Mu) && iscell(green) && ~iscell(Mv)
    I = ffm(Mu{1},X,green{1},k,Y,Mv,tol) + ...
        ffm(Mu{2},X,green{2},k,Y,Mv,tol) + ...
        ffm(Mu{3},X,green{3},k,Y,Mv,tol) ;
    
elseif ~iscell(Mu) && iscell(green) && iscell(Mv)
    I = ffm(Mu,X,green{1},k,Y,Mv{1},tol) + ...
        ffm(Mu,X,green{2},k,Y,Mv{2},tol) + ...
        ffm(Mu,X,green{3},k,Y,Mv{3},tol) ;
    
elseif iscell(Mu) && iscell(green) && iscell(Mv)
    I = ffm(Mu{1},X,green{2},k,Y,Mv{3},tol) - ...
        ffm(Mu{1},X,green{3},k,Y,Mv{2},tol) + ...
        ffm(Mu{2},X,green{3},k,Y,Mv{1},tol) - ...
        ffm(Mu{2},X,green{1},k,Y,Mv{3},tol) + ...
        ffm(Mu{3},X,green{1},k,Y,Mv{2},tol) - ...
        ffm(Mu{3},X,green{2},k,Y,Mv{1},tol) ;
else
    error('domIntegral7.m : unavailable case')
end
end



% function I = domIntegral7(data)
% %+========================================================================+
% %|                                                                        |
% %|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
% %|           openDom is part of the GYPSILAB toolbox for Matlab           |
% %|                                                                        |
% %| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
% %| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
% %| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
% %| LICENCE   : This program is free software, distributed in the hope that|
% %| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
% %| redistribute and/or modify it under the terms of the GNU General Public|
% %| License, as published by the Free Software Foundation (version 3 or    |
% %| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
% %| is available, please contact us to activate a "pay for remove" option. |
% %| CONTACT   : matthieu.aussal@polytechnique.edu                          |
% %|             francois.alouges@polytechnique.edu                         |
% %| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
% %|                                                                        |
% %| Please acknowledge the gypsilab toolbox in programs or publications in |
% %| which you use it.                                                      |
% %|________________________________________________________________________|
% %|   '&`   |                                                              |
% %|    #    |   FILE       : domIntegral7.m                                |
% %|    #    |   VERSION    : 0.40                                          |
% %|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
% %|  ( # )  |   CREATION   : 14.03.2017                                    |
% %|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
% %| ( === ) |   SYNOPSIS   : Numerical integation with 6 input arguments   |
% %|  `---'  |                                                              |
% %+========================================================================+
% 
% %%% BLOCK-MATRIX BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
% if (data{7} >= 1)
%     % Domain with quadrature
%     Xdom   = data{1};
%     [X,Wx] = Xdom.qud;
%     Nx     = size(X,1);
%     Wx     = spdiags(Wx,0,Nx,Nx);
%     
%     % Domain with quadrature
%     Ydom   = data{2};
%     [Y,Wy] = Ydom.qud;
%     Ny     = size(Y,1);
%     Wy     = spdiags(Wy,0,Ny,Ny);
%     
%     % Finite element matrix with integration
%     u  = data{3};
%     Mu = u.uqm(Xdom);
%     if iscell(Mu)
%         Mu{1} = Mu{1}' * Wx;
%         Mu{2} = Mu{2}' * Wx;
%         Mu{3} = Mu{3}' * Wx;
%     else
%         Mu = Mu' * Wx;
%     end
%     
%     % Green kernel
%     green = data{4};
%     
%     % Finite element matrix with integration
%     v  = data{5};
%     Mv = v.uqm(Ydom);
%     if iscell(Mv)
%         Mv{1} = Wy * Mv{1};
%         Mv{2} = Wy * Mv{2};
%         Mv{3} = Wy * Mv{3};
%     else
%         Mv = Wy * Mv;
%     end
%     
%     % Accuracy
%     tol = data{6};
%     
%     % Block
%     Nblk = data{7};
%     
%     % H-Matrix Integration
%     if iscell(Mu) && ~iscell(green) && ~iscell(Mv)
%         I{1} = hbm(u.unk,v.unk,Mu{1},X,green,Y,Mv,tol,Nblk);
%         I{2} = hbm(u.unk,v.unk,Mu{2},X,green,Y,Mv,tol,Nblk);
%         I{3} = hbm(u.unk,v.unk,Mu{3},X,green,Y,Mv,tol,Nblk);
%         
%     elseif ~iscell(Mu) && iscell(green) && ~iscell(Mv)
%         I{1} = hbm(u.unk,v.unk,Mu,X,green{1},Y,Mv,tol,Nblk);
%         I{2} = hbm(u.unk,v.unk,Mu,X,green{2},Y,Mv,tol,Nblk);
%         I{3} = hbm(u.unk,v.unk,Mu,X,green{3},Y,Mv,tol,Nblk);
%         
%     elseif ~iscell(Mu) && ~iscell(green) && iscell(Mv)
%         I{1} = hbm(u.unk,v.unk,Mu,X,green,Y,Mv{1},tol,Nblk);
%         I{2} = hbm(u.unk,v.unk,Mu,X,green,Y,Mv{2},tol,Nblk);
%         I{3} = hbm(u.unk,v.unk,Mu,X,green,Y,Mv{3},tol,Nblk);
%         
%     else
%         I = hbm(u.unk,v.unk,Mu,X,green,Y,Mv,tol,Nblk);
%     end
%     
%     
% %%% SCSD ELEMENT OPERATOR --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
% elseif (data{7} < 1)
%     % Domain with quadrature
%     Xdom   = data{1};
%     [X,Wx] = Xdom.qud;
%     Nx     = size(X,1);
%     Wx     = spdiags(Wx,0,Nx,Nx);
%     
%     % Domain with quadrature
%     Ydom   = data{2};
%     [Y,Wy] = Ydom.qud;
%     Ny     = size(Y,1);
%     Wy     = spdiags(Wy,0,Ny,Ny);
%     
%     % Finite element matrix with integration
%     u  = data{3};
%     Mu = u.uqm(Xdom);
%     if iscell(Mu)
%         Mu{1} = Mu{1}' * Wx;
%         Mu{2} = Mu{2}' * Wx;
%         Mu{3} = Mu{3}' * Wx;
%     else
%         Mu = Mu' * Wx;
%     end
%     
%     % Green kernel
%     green = data{4};
%     
%     % Wave number
%     k = data{5};
%     
%     % Finite element matrix with integration
%     v  = data{6};
%     Mv = v.uqm(Ydom);
%     if iscell(Mv)
%         Mv{1} = Wy * Mv{1};
%         Mv{2} = Wy * Mv{2};
%         Mv{3} = Wy * Mv{3};
%     else
%         Mv = Wy * Mv;
%     end
%     
%     % Accuracy
%     tol = data{7};
%     
%     % SCSD integration
%     if iscell(Mu) && ~iscell(green) && ~iscell(Mv)
%         I{1} = scm(Mu{1},X,green,k,Y,Mv,tol);
%         I{2} = scm(Mu{2},X,green,k,Y,Mv,tol);
%         I{3} = scm(Mu{3},X,green,k,Y,Mv,tol);
%         
%     elseif ~iscell(Mu) && iscell(green) && ~iscell(Mv)
%         I{1} = scm(Mu,X,green{1},k,Y,Mv,tol);
%         I{2} = scm(Mu,X,green{2},k,Y,Mv,tol);
%         I{3} = scm(Mu,X,green{3},k,Y,Mv,tol);
%         
%     elseif ~iscell(Mu) && ~iscell(green) && iscell(Mv)
%         I{1} = scm(Mu,X,green,k,Y,Mv{1},tol);
%         I{2} = scm(Mu,X,green,k,Y,Mv{2},tol);
%         I{3} = scm(Mu,X,green,k,Y,Mv{3},tol);
%         
%     elseif ~iscell(Mu) && iscell(green) && iscell(Mv)
%         I = scm(Mu,X,green{1},k,Y,Mv{1},tol) + ...
%             scm(Mu,X,green{2},k,Y,Mv{2},tol) + ...
%             scm(Mu,X,green{3},k,Y,Mv{3},tol) ;
%         
%     elseif iscell(Mu) && ~iscell(green) && iscell(Mv)
%         I = scm(Mu,X,green,k,Y,Mv,tol);
%         
%     elseif iscell(Mu) && iscell(green) && ~iscell(Mv)
%         I = scm(Mu{1},X,green{1},k,Y,Mv,tol) + ...
%             scm(Mu{2},X,green{2},k,Y,Mv,tol) + ...
%             scm(Mu{3},X,green{3},k,Y,Mv,tol) ;
%         
%     else
%         I = scm(Mu,X,green,k,Y,Mv,tol);
%     end
%     
%     
% else
%     error('domIntegral7.m : unavailable case')
% end
% end
