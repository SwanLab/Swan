function Mh = hmxBuilder(X,Y,M,tol)
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
%|    #    |   FILE       : hmxBuilder.m                                  |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Particles builder with low-rank approximation |
%|  `---'  |                for full, sparse and handle function          |
%+========================================================================+

% Initialisation
Mh = hmx(X,Y,tol);

% Admissibility
[isfar,Xdim,Ydim] = hmxFar(Mh);

% Compression for far distances
if ~issparse(M) && isfar
    % ACA with partial pivoting for handle function
    if isa(M,'function_handle')
        [A,B,flag] = hmxACA(X,Y,M,tol);
        
    % ACA with total pivoting for full matrix
    else
        [A,B,flag] = hmxACA(M,tol);
    end
    
% Compression for empty matrix
elseif issparse(M) && (nnz(M) == 0)
    A    = zeros(size(M,1),0);
    B    = zeros(0,size(M,2));
    flag = 1;
    
% No compression    
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
    
    % Handle function
    if isa(M,'function_handle')
        [I,J]  = ndgrid(1:size(X,1),1:size(Y,1));
        Mh.dat = M(X(I,:),Y(J,:));
        Mh.dat = reshape(Mh.dat,size(X,1),size(Y,1));
    
    % Full or sparse matrix
    else
        Mh.dat = M;
    end
    
    
%%% H-Matrix (recursion)
else
    % Type
    Mh.typ = 0; 
    
    % Subdivision for X
    [I1,I2] = hmxSubdivide(X,Xdim);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Y,Ydim);
    Mh.col  = {I1 , I2 , I1 , I2};
    
    % Single class
    if isa(X,'single')
        for i = 1:4
            Mh.row{i} = single(Mh.row{i});
            Mh.col{i} = single(Mh.col{i});
        end
    end

    % H-Matrix (recursion)
    for i = 1:4
        % Coordinates
        Xi = X(Mh.row{i},:);
        Yi = Y(Mh.col{i},:);
        
        % Partial pivoting
        if isa(M,'function_handle')
            Mh.chd{i} = hmxBuilder(Xi,Yi,M,tol);
            
        % Total pivoting    
        elseif isnumeric(M)
            Mi = M(Mh.row{i},Mh.col{i});
            Mh.chd{i} = hmxBuilder(Xi,Yi,Mi,tol);
            
        else
            error('hmxBuilder.m : unavailable case')
        end
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    
end
end
