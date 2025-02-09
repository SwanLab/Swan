function Ml = bmmCat(dim,Ml,Mr)
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
%|    #    |   FILE       : bmmCat.m                                      |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

%%% Vertical concatenation 
if (dim == 1) 
    % Security check
    if (norm(cell2mat(Ml.col)-cell2mat(Mr.col),'inf') > 0)
        error('bmm.m : unable to vertcat, please chek input dimension')
    end
   
    % Translate row indices for right matrix 
    for i = 1:length(Mr.row)
        Mr.row{i} = Mr.row{i} + size(Ml,1);
    end
    
    % Update left matrix
    Ml.row = [Ml.row;Mr.row];
    Ml.blk = [Ml.blk;Mr.blk];
    

%%% Horizontal concatenation 
elseif (dim == 2)
    % Security check
    if (norm(cell2mat(Ml.row)-cell2mat(Mr.row),'inf') > 0)
        error('bmm.m : unable to horzcat, please chek input dimension')
    end
    
    % Translate column indices for right matrix 
    for j = 1:length(Mr.col)
        Mr.col{j} = Mr.col{j} + size(Ml,2);
    end
    
    % Update left matrix
    Ml.col = [Ml.col;Mr.col];
    Ml.blk = [Ml.blk,Mr.blk];
end
end
