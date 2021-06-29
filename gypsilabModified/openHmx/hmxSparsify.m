function Ms = hmxSparsify(Mh,Ms)
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
%|    #    |   FILE       : hmxSparsify.m                                 |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix to sparse matrix             |
%|  `---'  |                                                              |
%+========================================================================+
    
%%% H-Matrix (recursion)
if (Mh.typ == 0)
    for i = 1:4
        Ms(Mh.row{i},Mh.col{i}) = hmxSparsify(Mh.chd{i},Ms(Mh.row{i},Mh.col{i}));
    end
    
%%% Compressed leaf
elseif (Mh.typ == 1)
    [I,J] = find(Ms);
    A     = Mh.dat{1}(I,:);
    B     = Mh.dat{2}(:,J);
    Ms(sub2ind(size(Ms),I,J)) = sum(A.*B.',2);    
    
%%% Full leaf
elseif (Mh.typ == 2)
    I     = find(Ms);
    Ms(I) = Mh.dat(I);

%%% Unknown type
else
    error('hmxSparsify.m : unavailable case')
end
end
