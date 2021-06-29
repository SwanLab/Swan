function Mh = hmxInv(Mh)
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
%|    #    |   FILE       : hmxInv.m                                      |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Inversion of H-Matrix based on Schur          |
%|  `---'  |                complement                                    |
%+========================================================================+

% Check invertibility
if (size(Mh,1) ~= size(Mh,2)) || ~isequal(Mh.pos{1},Mh.pos{2})
    error('hmxInv.m : matrix must be invertible.')
end

%%% H-Matrix (bloc recursion)
if (Mh.typ == 0)
    % Am1 -> M11 
    Mh.chd{1} = hmxInv(Mh.chd{1});
    
    % - Am1 * B -> X12
    X12 = - Mh.chd{1} * Mh.chd{2};
    
    % C * Am1 -> X21
    X21 = Mh.chd{3} * Mh.chd{1};
    
    % D - C * Am1 * B = S -> M22  (schur complement)
    Mh.chd{4} = Mh.chd{4} + Mh.chd{3} * X12;
    
    % S^-1 -> M22
    Mh.chd{4} = hmxInv(Mh.chd{4});
    
    % - Am1 * B * S -> M12
    Mh.chd{2} = X12 * Mh.chd{4};
    
    % Am1 + Am1 * B * S * C * Am1 -> M11
    Mh.chd{1} = Mh.chd{1} - Mh.chd{2} * X21;
    
    % - S * C * Am1 -> M21
    Mh.chd{3} = - Mh.chd{4} * X21;
    
    % Fusion
    Mh = hmxFusion(Mh);
    
%%% Compressed leaf    
elseif (Mh.typ == 1)
    error('hmxInv : unavailable case')
    
%%% Full leaf    
elseif (Mh.typ == 2)
    Mh.dat = inv(Mh.dat);

%%% Unknown type
else
    error('hmxInv.m : unavailable case')
end
end
