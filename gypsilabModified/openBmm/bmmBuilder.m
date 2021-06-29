function Mb = bmmBuilder(data)
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
%|    #    |   FILE       : bmmBuilder.m                                  |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

%%% Block matrix
if isa(data,'bmm')
    Mb = data;
  
    
%%% Cells
elseif iscell(data)
    % Dimension for each block
    I = zeros(size(data));
    J = zeros(size(data));
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            I(i,j) = size(data{i,j},1);
            J(i,j) = size(data{i,j},2);
        end
    end
    
    % Security
    if (norm(diff(I,[],2),'inf') > 0)
        error('bmmBuilder.m : unavailable dimensions, please check input')
    end
    if (norm(diff(J,[],1),'inf') > 0)
        error('bmmBuilder.m : unavailable dimensions, please check input')
    end
    
    % Loop for row
    Mb     = bmm();
    Mb.row = mat2cell((1:sum(I(:,1)))',I(:,1),1);
    Mb.col = mat2cell((1:sum(J(1,:)))',J(1,:)',1);
    Mb.blk = data;
    
    
%%% Single block
else
    Mb     = bmm();
    Mb.row = {(1:size(data,1))'};
    Mb.col = {(1:size(data,2))'};
    Mb.blk = {data};
end
end
