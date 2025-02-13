function Mh = hmxDouble(Mh)
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
%|    #    |   FILE       : hmxDouble.m                                   |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix to double precision          |
%|  `---'  |                                                              |
%+========================================================================+

% Position
Mh.pos = {double(Mh.pos{1}),double(Mh.pos{2})};

%%% H-Matrix (recursion)
if (Mh.typ == 0)
    for i = 1:4
        Mh.chd{i} = hmxDouble(Mh.chd{i});
        Mh.row{i} = double(Mh.row{i});
        Mh.col{i} = double(Mh.col{i});
    end
    Mh = hmxFusion(Mh);
    
%%% Compressed leaf
elseif (Mh.typ == 1)
    Mh.dat{1} = double(Mh.dat{1});
    Mh.dat{2} = double(Mh.dat{2});
    
%%% Full leaf
elseif (Mh.typ == 2)
    Mh.dat = double(Mh.dat);
    
%%% Unknown type
else
    error('hmxDouble.m : unavailable case')
end
end
