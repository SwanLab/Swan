function Mh = hmxBuilderFem(Xunk,Yunk,Mx,X,green,Y,My,tol)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : hmxBuilderFem.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Finite element builder with low-rank          |
%|  `---'  |                approximation for handle function             |
%+========================================================================+

% Initialisation
Mh = hmx(Xunk,Yunk,tol);

% Admissibility
[isfar,Xdim,Ydim] = hmxFar(Mh);

% Compression for far distances
if isfar
    % ACA compression
    if iscell(green)
        A = cell(1,3);
        B = cell(1,3);
        [A{1},B{1},flag1] = hmxACA(X,Y,green{1},tol);
        [A{2},B{2},flag2] = hmxACA(X,Y,green{2},tol);
        [A{3},B{3},flag3] = hmxACA(X,Y,green{3},tol);
        flag = flag1 * flag2 * flag3;
    else
        [A,B,flag] = hmxACA(X,Y,green,tol);
    end
    
    % Update
    if flag
        % Summation
        if iscell(Mx) && iscell(A) && iscell(My)
            A = [ Mx{1}*A{2} , - Mx{1}*A{3} , ...
                Mx{2}*A{3} , - Mx{2}*A{1} , ...
                Mx{3}*A{1} , - Mx{3}*A{2} ] ;
            B = [ B{2}*My{3} ; B{3}*My{2} ; ...
                B{3}*My{1} ; B{1}*My{3} ; ...
                B{1}*My{2} ; B{2}*My{1} ] ;
            
        elseif iscell(Mx) && iscell(A) && ~iscell(My)
            A = [Mx{1}*A{1} , Mx{2}*A{2} , Mx{3}*A{3}];
            B = [B{1}*My    ; B{2}*My    ; B{3}*My ];
            
        elseif iscell(Mx) && ~iscell(A) && iscell(My)
            A = [Mx{1}*A , Mx{2}*A , Mx{3}*A];
            B = [B*My{1} ; B*My{2} ; B*My{3} ];
            
        elseif ~iscell(Mx) && iscell(A) && iscell(My)
            A = [Mx*A{1}    , Mx*A{2}    , Mx*A{3}];
            B = [B{1}*My{1} ; B{2}*My{2} ; B{3}*My{3}];
            
        else
            A = Mx * A;
            B = B * My;
        end
        
        % Recompression
        [A,B] = hmxQRSVD(A,B,tol);
    end
    
else
    flag = 0;
end

    
%%% Compression
if flag
    % Type
    Mh.typ = 1;
    
    % Low-rank
    Mh.dat = {A,B};


%%%% Full or sparse for smallest box (stopping criterion)
elseif (min(size(Mh)) < 100)
    % Type
    Mh.typ = 2;

    % Quadrature matrix
    Nx    = size(X,1);
    Ny    = size(Y,1);
    [I,J] = ndgrid(1:Nx,1:Ny);
    if iscell(green)
        Gxy{1} = reshape(green{1}(X(I,:),Y(J,:)),Nx,Ny);
        Gxy{2} = reshape(green{2}(X(I,:),Y(J,:)),Nx,Ny);
        Gxy{3} = reshape(green{3}(X(I,:),Y(J,:)),Nx,Ny);
    else
        Gxy = reshape(green(X(I,:),Y(J,:)),Nx,Ny);
    end
    
    % Matrix integration
    Mh.dat = femMultiplyCell(Mx,Gxy,My);
    

%%% H-Matrix (recursion)
else
    % Type
    Mh.typ = 0;
    
    % Subdivision for X
    [I1,I2] = hmxSubdivide(Xunk,Xdim);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Yunk,Ydim);
    Mh.col  = {I1 , I2 , I1 , I2};

    % H-Matrix (recursion)
    for i = 1:4
        % Dof indices
        Ir = Mh.row{i};
        Ic = Mh.col{i};
        
        % Fem matrix subdivision and quadratures points indices
        [Mxchd,Ix] = femSubdivideCell(Mx,Ir,'left');
        [Mychd,Iy] = femSubdivideCell(My,Ic,'right');

        % Recursion
        Mh.chd{i} = hmxBuilderFem(Xunk(Ir,:),Yunk(Ic,:),...
            Mxchd,X(Ix,:),green,Y(Iy,:),Mychd,tol);
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    
end
end
