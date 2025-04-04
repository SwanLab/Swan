function D = hmxDiag(Mh,I,J)
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
%|    #    |   FILE       : hmxDiag.m                                     |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Extract diagonal from H-Matrix                |
%|  `---'  |                                                              |
%+========================================================================+
    
%%% H-Matrix (recursion)
if (Mh.typ == 0)
    D = sparse(size(Mh,1),size(Mh,2));
    for i = 1:4
        D(Mh.row{i},Mh.col{i}) = hmxDiag(Mh.chd{i},I(Mh.row{i}),J(Mh.col{i}));
    end
    
%%% Compressed leaf
elseif (Mh.typ == 1)
    % Diagonal indices
    [ia,ib] = ismember(I,J);
    idx     = find(ia);
    jdx     = ib(ia);

    % Diagonal vector
    D = sum(Mh.dat{1}(idx,:) .* (Mh.dat{2}(:,jdx)).',2);
    
    % Diagonal sparse matrix
    D = sparse(idx,jdx,D,size(Mh,1),size(Mh,2));
        
%%% Full leaf
elseif (Mh.typ == 2)
    % Diagonal indices
    [ia,ib] = ismember(I,J);
    idx     = find(ia);
    jdx     = ib(ia);
    
    % Diagonal vector
    D = Mh.dat(sub2ind(size(Mh.dat),idx,jdx));
    
    % Diagonal sparse matrix
    D = sparse(idx,jdx,D,size(Mh,1),size(Mh,2));
    
%%% Unknown type
else
    error('hmxDiag.m : unavailable case')
end
