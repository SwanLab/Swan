function Mh = hmxCtranspose(Mh)
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
%|    #    |   FILE       : hmxCtranspose.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Conjugate transposition of H-Matrix           |
%|  `---'  |                                                              |
%+========================================================================+

% Dimensions
Mh.pos = {Mh.pos{2} Mh.pos{1}};

%%% H-Matrix (recursion)
if (Mh.typ == 0)
    % Data
    tmp       = Mh.chd{2};
    Mh.chd{1} = hmxCtranspose(Mh.chd{1});    
    Mh.chd{2} = hmxCtranspose(Mh.chd{3});
    Mh.chd{3} = hmxCtranspose(tmp);
    Mh.chd{4} = hmxCtranspose(Mh.chd{4});

    % Indices
    I      = [1 3 2 4];
    tmp    = Mh.row;
    Mh.row = Mh.col(I);
    Mh.col = tmp(I);
    
%%% Compressed leaf
elseif (Mh.typ == 1)
    Mh.dat = {Mh.dat{2}' , Mh.dat{1}'};
    
%%% Full leaf
elseif (Mh.typ == 2)
    Mh.dat = Mh.dat';

%%% Unknown type
else
    error('hmxCtranspose.m : unavailable case')
end
end
