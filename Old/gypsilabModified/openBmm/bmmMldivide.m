function B = bmmMldivide(Mb,B)
%+========================================================================+
%|                                                                        |
%|               OPENBMM - LIBRARY FOR BLOCK MATRIX ALGEBRA               |
%|           openBmm is part of the GYPSILAB toolbox for Matlab           |
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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : bmmMldivide.m                                 |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Check block
n = size(Mb.blk,1);
if (size(Mb.blk,2) ~= n)
    error('bmmMldivide.m : matrix blocks must agree.')
end

% Check indices
if (norm(cell2mat(Mb.row)-cell2mat(Mb.col),'inf') > 0)
    error('bmmMldivide.m : matrix indices must agree.')
end

% Is triangular Block-Matrix
if (n > 1)
    % Initialize test
    islower = true;
    isupper = true;
    
    % Double loop for extra-diagonal with symetry
    for i = 1:n
        for j = i+1:n
            % Lower
            M = Mb.blk{i,j};
            if isa(M,'hmx')
                islower = islower && (M.typ == 1) && (isempty(M.dat{1}));
            else
                islower = islower && (nnz(M)==0);
            end
            
            % Upper
            M = Mb.blk{j,i};
            if isa(M,'hmx')
                isupper = isupper && (M.typ == 1) && (isempty(M.dat{1}));
            else
                isupper = isupper && (nnz(M)==0);
            end
        end
    end

% Is mono block matrix   
else
    if isa(B,'bmm')
        B.blk{1} = Mb.blk{1}\B.blk{1};
    else
        [L,U] = lu(Mb.blk{1});
        B(Mb.col{1},:) = U \ (L \ B(Mb.row{1},:));
    end
    return
end    
    

%%% System is islower triangular and RHS is numeric
if islower && isnumeric(B)
    % For each row
    for k = 1:n
        % Load Lkk block
        Lkk = Mb.blk{k,k};
        
        % Load Bk rhs
        Bk = B(Mb.row{k},:);
        
        % Bk = Bk - sum_{l=1}^{k-1} Lkl Bl
        for l = 1:k-1
            Bk  = Bk - Mb.blk{k,l} * B(Mb.row{l},:);
        end
        
        % Solve linear system 
        B(Mb.col{k},:) = Lkk \ Bk;
    end

    
%%% System is isupper triangular and RHS is numeric
elseif isupper && isnumeric(B)
    % For each row
    for k = (n:-1:1)
        % Load Ukk block
        Ukk = Mb.blk{k,k};
        
        % Load Bk rhs
        Bk = B(Mb.row{k},:);
        
        % Bk = Bk - sum_{l=k+1}^{n} Ukl Bl
        for l = k+1:n
            Bk  = Bk - Mb.blk{k,l} * B(Mb.row{l},:);
        end
        
        % Solve linear system
        B(Mb.col{k},:) = Ukk \ Bk;
    end

    
%%% System is islower triangular and RHS is block
elseif islower && isa(B,'bmm')
    % For each lhs row
    for k = 1:n
        % Load Lkk block
        Lkk = Mb.blk{k,k};
        
        % For each rhs column
        for l = 1:size(B.blk,2)
            % Load Bkl rhs
            Bkl = B.blk{k,l};
            
            % Bkl = Bkl - sum_{m=1}^{k-1} Lkm Bml
            for m = 1:k-1
                Bkl = Bkl - Mb.blk{k,m} * B.blk{m,l};
            end
            
            % Solve linear system
            B.blk{k,l} = Lkk \ Bkl;
        end
    end

    
%%% System is isupper triangular and RHS is block
elseif isupper && isa(B,'bmm')
    % For each lhs row
    for k = (n:-1:1)
        % Load Ukk block
        Ukk = Mb.blk{k,k};
        
        % For each rhs column
        for l = 1:size(B.blk,2)
            % Load Bkl rhs
            Bkl = B.blk{k,l};
            
            % Bkl = Bkl - sum_{m=k+1}^{n} Ukm Bml
            for m = k+1:length(Mb.col)
                Bkl = Bkl - Mb.blk{k,m} * B.blk{m,l};
            end
            
            % Solve linear system
            B.blk{k,l} = Ukk \ Bkl;
        end
    end
    
    
%%% System is not triangular
else
    % LU factorization
    [Lb,Ub] = lu(Mb);
    
    % Solve
    B = Ub \ (Lb \ B);
end        
end
