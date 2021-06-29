function [Mb,Ub] = bmmLU(Mb)
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
%|    #    |   FILE       : bmmLU.m                                       |
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
    error('bmmLU.m : matrix blocks must agree.')
end

% Check indices
if (norm(cell2mat(Mb.row)-cell2mat(Mb.col),'inf') > 0)
    error('bmmLU.m : matrix indices must agree.')
end

% Initialization
Ub     = bmm();
Ub.row = Mb.row;
Ub.col = Mb.col;

% For diagonal dimension
for k = 1:n
    %%%%%%%%%%%%%%%%%%%%%%%%%%% DIAGONAL BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Mkk block
    Mkk = Mb.blk{k,k};
    
    % Mkk = Mkk - sum_{l=1}^{k-1} Lkl Ulk 
    for l = 1:k-1
        Mkk = Mkk - Mb.blk{k,l} * Ub.blk{l,k};
    end
    
    % LU factorization : Lkk*Ukk = Mkk - sum_{l=1}^{k-1} Lkl Ulk 
    [Mkk,Ukk] = lu(Mkk);
    
    % Save Lkk and Ukk block
    Mb.blk{k,k} = Mkk;
    Ub.blk{k,k} = Ukk;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROW BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop on row blocks
    for i = k+1:n
        % Load Mik block
        Mik = Mb.blk{i,k};
        
        % Mik = Mik - sum_{l=1}^{k-1} Lil Ulk
        for l = 1:k-1
            Mik = Mik - Mb.blk{i,l} * Ub.blk{l,k};
        end
        
        % Solve upper system : Lik*Ukk = Mik - sum_{l=1}^{k-1} Lil Ulk
        Mb.blk{i,k} = Mik / Ukk;
        
        % Save Uik block (empty)
        dim = size(Mb.blk{i,k});
        Ub.blk{i,k} = sparse(dim(1),dim(2));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% COLUMN BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop on column blocks
    for j = k+1:n
        % Load Mkj block
        Mkj = Mb.blk{k,j};
        
        % Mkj = Mkj - sum_{l=1}^{k-1} Lkl Ulj
        for l = 1:k-1
            Mkj = Mkj - Mb.blk{k,l} * Ub.blk{l,j};
        end
        
        % Solve lower system : Lkk*Ukj = Mkj - sum_{l=1}^{k-1} Llk Ulj
        Ub.blk{k,j} = Mkk \ Mkj;
        
        % Save Lkj block (empty)
        dim = size(Ub.blk{k,j});
        Mb.blk{k,j} = sparse(dim(1),dim(2));
    end
        
end
end
