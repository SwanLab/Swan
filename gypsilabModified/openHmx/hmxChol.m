function Mh = hmxChol(Mh)
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
%|    #    |   FILE       : hmxChol.m                                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Cholesky factorization of H-Matrix            |
%|  `---'  |                                                              |
%+========================================================================+

% Check invertibility
if (size(Mh,1) ~= size(Mh,2)) || ~isequal(Mh.pos{1},Mh.pos{2})
    error('hmxChol.m : matrix must be invertible.')
end

%%% H-Matrix (recursion)
if (Mh.typ == 0)
    % U11 -> M11
    Mh.chd{1} = hmxChol(Mh.chd{1});

    % U12 -> U11' \ M12
    Mh.chd{2} = Mh.chd{1}'\Mh.chd{2};
    
    % U21 -> 0
    Mh.chd{3} = zeros(Mh.chd{3});
    
    % M22 -> M22 - U12'*U12
    Mh.chd{4} = Mh.chd{4} - Mh.chd{2}' * Mh.chd{2};
%     Mh.chd{4} = plusmtimes(Mh.chd{4},-1,Mh.chd{2}',Mh.chd{2});
    
    % U22 -> M22
    Mh.chd{4} = hmxChol(Mh.chd{4});
    
    % Fusion
    Mh = hmxFusion(Mh);

%%% Compressed leaf
elseif (Mh.typ == 1)
    error('hmxChol : unavailable case')
    
%%% Full leaf
elseif (Mh.typ == 2)
    Mh.dat = chol(Mh.dat);    
   
%%% Unknown type
else
    error('hmxChol.m : unavailable case')
end
end
