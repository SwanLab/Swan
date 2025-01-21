function M = bmmMtimes(Ml,Mr)
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
%|    #    |   FILE       : bmmMtimes.m                                   |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Check dimensions
if (size(Ml,2) ~= size(Mr,1))
    error('bmmMtimes.m : matrix dimensions must agree.')
end


%%% Block-Matrix * Block-Matrix
if isa(Ml,'bmm') && isa(Mr,'bmm')
    % Check block 
    if (size(Ml.blk,2) ~= size(Mr.blk,1))
        error('bmmMtimes.m : matrix blocks must agree.')
    end
    
    % Check indices
    if (norm(cell2mat(Ml.col)-cell2mat(Mr.row),'inf') > 0)
        error('bmmMtimes.m : matrix indices must agree.')
    end

    % Initialize
    M     = bmm();
    M.row = Ml.row;
    M.col = Mr.col;
    M.blk = cell(length(Ml.row),length(Mr.col));
        
    % Double loop
    for i = 1:length(Ml.row)
        for j = 1:length(Mr.col)
            % Initialize block
            M.blk{i,j} = sparse(length(M.row{i}),length(M.col{j}));
                        
            % Block product
            for k = 1:length(Ml.col)
                M.blk{i,j} = M.blk{i,j} + Ml.blk{i,k}*Mr.blk{k,j};
            end
        end
    end
    
    
%%% Block-Matrix * Full --> Full
elseif isa(Ml,'bmm') && isnumeric(Mr)
    % Initialization
    M = zeros(size(Ml,1),size(Mr,2));

    % For each block
    for i = 1:length(Ml.row)
        for j = 1:length(Ml.col)
            M(Ml.row{i},:) = M(Ml.row{i},:) + Ml.blk{i,j}*Mr(Ml.col{j},:);
        end
    end
    
    
%%% Full * Block-Matrix --> Full
elseif isnumeric(Ml) && isa(Mr,'bmm')
    % Initialization
    M = zeros(size(Ml,1),size(Mr,2));

    % For each block
    for i = 1:length(Mr.row)
        for j = 1:length(Mr.col)
            M(:,Mr.col{j}) = M(:,Mr.col{j}) + Ml(:,Mr.row{i})*Mr.blk{i,j};
        end
    end
    
    
%%% Unavailable    
else
    error('bmmMtimes.m : unavailable case')
end
end
