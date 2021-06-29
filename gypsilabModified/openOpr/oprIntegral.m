function Mh = oprIntegral(opr,k,Xdom,Xfem,Ydom,Yfem,tol)
%+========================================================================+
%|                                                                        |
%|            OPENOPR - LIBRARY FOR SPECIFIC OPERATORS IN BEM             |
%|           openOpr is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : oprIntegral.m                                 |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 8 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

% Domain X with quadrature
[X,Wx] = Xdom.qud;
Xnrm   = Xdom.qudNrm;
Nx     = size(X,1);
Wx     = spdiags(Wx,0,Nx,Nx);

% Domain Y with quadrature
[Y,Wy] = Ydom.qud;
Ynrm   = Ydom.qudNrm;
Ny     = size(Y,1);
Wy     = spdiags(Wy,0,Ny,Ny);


%%%%%%%%%% HELMHOLTZ : Single layer (S)
if strcmp(opr,'S') 
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = (1/(4*pi)) .* Xfem.uqm(Xdom)' * Wx;
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = Wy * Yfem.uqm(Ydom);
    
    
%%%%%%%%%% HELMHOLTZ : double layer (D)  
elseif strcmp(opr,'D')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = (1/(4*pi)) .* Xfem.uqm(Xdom)' * Wx;
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = Wy * Yfem.uqm(Ydom);
    

%%%%%%%%%% HELMHOLTZ : double layer transpose (Dt)    
elseif strcmp(opr,'Dt')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = (1/(4*pi)) .* Xfem.uqm(Xdom)' * Wx;
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = Wy * Yfem.uqm(Ydom);
 
    
%%%%%%%%%% HELMHOLTZ : Hypersingular (H)
elseif strcmp(opr,'H')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = [uqm(ntimes(Xfem),Xdom),uqm(nxgrad(Xfem),Xdom)];
    cte  = 1/(4*pi) .* [k^2,k^2,k^2,-1,-1,-1];
    for n = 1:length(Mx)
        Mx{n} = cte(n) * Mx{n}' * Wx;
    end
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = [uqm(ntimes(Yfem),Ydom),uqm(nxgrad(Yfem),Ydom)];
    for n = 1:length(My)
        My{n} = Wy * My{n};
    end
    

%%%%%%%%%% MAXWELL : EFIE (T)
elseif strcmp(opr,'T')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = [uqm(Xfem,Xdom),{uqm(div(Xfem),Xdom)}];
    cte  = 1/(4*pi) .* [1i*k,1i*k,1i*k,-1i/k];
    for n = 1:length(Mx)
        Mx{n} = cte(n) * Mx{n}' * Wx;
    end
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = [uqm(Yfem,Ydom),{uqm(div(Yfem),Ydom)}];
    for n = 1:length(My)
        My{n} = Wy * My{n};
    end

  
%%%%%%%%%% MAXWELL : MFIE (nxK)
elseif strcmp(opr,'nxK')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = uqm(nx(Xfem),Xdom);
    cte  = 1/(4*pi);
    for n = 1:length(Mx)
        Mx{n} = cte * Mx{n}' * Wx;
    end
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = uqm(Yfem,Ydom);
    for n = 1:length(My)
        My{n} = Wy * My{n};
    end
    

%%%%%%%%%% STOKES : Stokeslet (G)    
elseif strcmp(opr,'G')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = (1/(8*pi)) .* Xfem.uqm(Xdom)' * Wx;
        
    % Vectorial X 
    Zx   = sparse(size(Xunk,1),size(X,1));
    Mx   = [Mx Zx Zx ; Zx Mx Zx ; Zx Zx Mx];
    X    = [X ones(Nx,1) ; X 2*ones(Nx,1) ; X 3*ones(Nx,1)];
    Xnrm = [Xnrm ; Xnrm ; Xnrm];
    Xunk = [Xunk ; Xunk ; Xunk];
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = Wy * Yfem.uqm(Ydom); 

    % Vectorial Y
    Zy   = sparse(size(Y,1),size(Yunk,1));
    My   = [My Zy Zy ; Zy My Zy ; Zy Zy My];
    Y    = [Y ones(Ny,1) ; Y 2*ones(Ny,1) ; Y 3*ones(Ny,1)];
    Ynrm = [Ynrm ; Ynrm ; Ynrm];
    Yunk = [Yunk ; Yunk ; Yunk];
  
       
%%%%%%%%%% STOKES : Stresslet (T)    
elseif strcmp(opr,'Ts')
    % Finite element X
    Xunk = Xfem.unk;
    Mx   = (-6/(8*pi)) .* Xfem.uqm(Xdom)' * Wx;
        
    % Vectorial X 
    Zx   = sparse(size(Xunk,1),size(X,1));
    Mx   = [Mx Zx Zx ; Zx Mx Zx ; Zx Zx Mx];
    X    = [X ones(Nx,1) ; X 2*ones(Nx,1) ; X 3*ones(Nx,1)];
    Xnrm = [Xnrm ; Xnrm ; Xnrm];
    Xunk = [Xunk ; Xunk ; Xunk];
    
    % Finite element Y
    Yunk = Yfem.unk;
    My   = Wy * Yfem.uqm(Ydom); 

    % Vectorial Y
    Zy   = sparse(size(Y,1),size(Yunk,1));
    My   = [My Zy Zy ; Zy My Zy ; Zy Zy My];
    Y    = [Y ones(Ny,1) ; Y 2*ones(Ny,1) ; Y 3*ones(Ny,1)];
    Ynrm = [Ynrm ; Ynrm ; Ynrm];
    Yunk = [Yunk ; Yunk ; Yunk];

    
%%%%%%%%%% Unknown    
else    
    error('oprIntegral.m : unknown operator')
end


%%%%%%%%%% Build operator
Mh = oprBuilderHmx(opr,k,Xunk,Mx,X,Xnrm,Yunk,My,Y,Ynrm,tol);    
end




%     % Y integration weigth vectorial
%     Wy = spdiags(Wy);
%     Wy = [Wy;Wy;Wy];
% 
%     % Builder
%     [Mh,Tun] = oprBuilderStk('Ts',Xunk,Mx,X,Xnrm,Yunk,My,Y,Ynrm,Wy,tol);
% 
%     % X Integration matrix 
%     Mx = (-6/(8*pi)) .* Xfem.uqm(Xdom)' * Wx;
%     Mx = [Mx Mx Mx ; Mx Mx Mx ; Mx Mx Mx];
%     
%     % Y integration matrix
%     My = Yfem.uqm(Ydom); 
%     My = [My My My ; My My My ; My My My];
%     
%     % Add correction
%     Mh = Mh - Mx * spdiags(Tun,0,3*Nx,3*Nx) * My;
