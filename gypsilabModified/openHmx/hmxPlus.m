function Ml = hmxPlus(Ml,Mr)
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
%|    #    |   FILE       : hmxPlus.m                                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Sum of H-Matrix with the rule                 |
%|  `---'  |                Full > H-Matrix > Compr                       |
%+========================================================================+

% Check dimensions
if sum(size(Ml) ~= size(Mr))
    error('hmxPlus.m : matrix dimensions must agree.')
end

%%% H-Matrix + H-Matrix -> H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    % H-Matrix + H-Matrix --> H-Matrix (recursion)
    if (Ml.typ == 0) && (Mr.typ == 0)
        % Construction
        for i = 1:4
            Ml.chd{i} = hmxPlus(Ml.chd{i},Mr.chd{i});
        end

        % Fusion
        Ml = hmxFusion(Ml);
    
    % H-Matrix + Compr -> H-Matrix    
    elseif (Ml.typ == 0) && (Mr.typ == 1)
        Ml = hmxPlusAB(Ml,Mr.dat{1},Mr.dat{2});
             
    % H-Matrix + Full -> Unavailable
    elseif (Ml.typ == 0) && (Mr.typ == 2)
        error('hmxPlus : unvailable case')

        
    % Compr + --- -> ---
    elseif (Ml.typ == 1)
        Ml = hmxPlusAB(Mr,Ml.dat{1},Ml.dat{2});
        
        
    % Full + H-Matrix -> Unavailable
    elseif (Ml.typ == 2) && (Mr.typ == 0)
        error('hmxPlus : unvailable case')
        
    % Full + Compr -> Full
    elseif (Ml.typ == 2) && (Mr.typ == 1)
        Ml = hmxPlusAB(Ml,Mr.dat{1},Mr.dat{2});

    % Full + Full -> Full                              
    elseif (Ml.typ == 2) && (Mr.typ == 2)
        Ml.dat = Ml.dat + Mr.dat;
        
        
    else
        error('hmxPlus : unvailable case')
    end

    
    
%%% H-Matrix + Matrix -> H-Matrix
elseif isa(Ml,'hmx') 
    Ml = Ml + hmx(Ml.pos{1},Ml.pos{2},Mr,Ml.tol);

    
%%% Matrix + H-Matrix -> Matrix
elseif isa(Mr,'hmx')
    Ml = hmx(Mr.pos{1},Mr.pos{2},Ml,Mr.tol) + Mr;

    
%%% Unavailable  
else
    error('hmxPlus.m : unavailable case')
end
end
