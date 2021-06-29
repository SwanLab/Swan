function Mh = hmxPlusMtimes(Mh,alpha,Ml,Mr)
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
%|    #    |   FILE       : hmxPlusMtimes.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Sum and product of H-Matrix with the rule     |
%|  `---'  |                                                              |
%+========================================================================+

% Check dimensions
if (size(Mh,1) ~= size(Ml,1))
    error('hmxPlusMtimes.m : matrix dimensions must agree.')
end
if (size(Mh,2) ~= size(Mr,2))
    error('hmxPlusMtimes.m : matrix dimensions must agree.')
end
if (size(Ml,2) ~= size(Mr,1))
    error('hmxPlusMtimes.m : matrix dimensions must agree.')
end

%%% H-Matrix + alpha * H-Matrix * H-Matrix --> H-Matrix
if isa(Mh,'hmx') && isa(Ml,'hmx') && isa(Mr,'hmx')
    
    % H-Matrix + H-Matrix * H-Matrix --> H-Matrix   (recursion)
    if (Mh.typ==0) && (Ml.typ==0) && (Mr.typ==0)
        % First bloc : M1 -> M1 + alpha * (Ml1 * Mr1 + Ml2 * Mr3)
        Mh.chd{1} = hmxPlusMtimes(Mh.chd{1},alpha,Ml.chd{1},Mr.chd{1});
        Mh.chd{1} = hmxPlusMtimes(Mh.chd{1},alpha,Ml.chd{2},Mr.chd{3});
       
        % Second bloc : M2 -> M2 + alpha * (Ml1 * Mr2 + Ml2 * Mr4) 
        Mh.chd{2} = hmxPlusMtimes(Mh.chd{2},alpha,Ml.chd{1},Mr.chd{2});
        Mh.chd{2} = hmxPlusMtimes(Mh.chd{2},alpha,Ml.chd{2},Mr.chd{4});
        
        % Third bloc : M3 -> M3 + alpha * (Ml3 * Mr1 + Ml4 * Mr3)
        Mh.chd{3} = hmxPlusMtimes(Mh.chd{3},alpha,Ml.chd{3},Mr.chd{1});
        Mh.chd{3} = hmxPlusMtimes(Mh.chd{3},alpha,Ml.chd{4},Mr.chd{3});
        
        % Fourth bloc : M4 -> M4 + alpha * (Ml3 * Mr2  + Ml4 * Mr4)
        Mh.chd{4} = hmxPlusMtimes(Mh.chd{4},alpha,Ml.chd{3},Mr.chd{2});
        Mh.chd{4} = hmxPlusMtimes(Mh.chd{4},alpha,Ml.chd{4},Mr.chd{4});
        
        % Fusion
        Mh = hmxFusion(Mh);
        
    % H-Matrix + H-Matrix * Compr --> H-Matrix
    elseif (Mh.typ==0) && (Ml.typ==0) && (Mr.typ==1)
        A  = alpha * ( Ml * Mr.dat{1} );
        B  = Mr.dat{2};
        Mh = hmxPlusAB(Mh,A,B);

    % H-Matrix + H-Matrix * Full --> unavailable
    elseif (Mh.typ==0) && (Ml.typ==0) && (Mr.typ==2)
        error('hmxPlusMtimes : unvailable case')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    % H-Matrix + Compr * H-Matrix --> H-Matrix
    elseif (Mh.typ==0) && (Ml.typ==1) && (Mr.typ==0)
        A  = alpha * Ml.dat{1}; 
        B  = Ml.dat{2} * Mr;
        Mh = hmxPlusAB(Mh,A,B);
       
    % H-Matrix + Compr * Compr --> H-Matrix
    elseif (Mh.typ==0) && (Ml.typ==1) && (Mr.typ==1) 
        A  = alpha * Ml.dat{1};
        B  = (Ml.dat{2} * Mr.dat{1}) * Mr.dat{2};
        Mh = hmxPlusAB(Mh,A,B);
        
    % H-Matrix + Compr * Full --> H-Matrix
    elseif (Mh.typ==0) && (Ml.typ==1) && (Mr.typ==2)
        error('hmxPlusMtimes : unvailable case')
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    % H-Matrix + Full * H-Matrix --> Full
    elseif (Mh.typ==0) && (Ml.typ==2)
        error('hmxPlusMtimes : unvailable case')
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
    % Compr + H-Matrix * H-Matrix --> H-Matrix
    elseif (Mh.typ==1) && (Ml.typ==0) && (Mr.typ==0)
        Mh.typ = 0;
        Mh.row = Ml.row;
        Mh.col = Mr.col;
        for i = 1:4
            Mh.chd{i} = hmx(Ml.chd{i}.pos{1},Mr.chd{i}.pos{2},...
                Mh.dat{1}(Ml.row{i},:) , Mh.dat{2}(:,Mr.col{i}) , ...
                Ml.tol ) ;
        end
        Mh.dat = [];
        Mh     = hmxPlusMtimes(Mh,alpha,Ml,Mr);

    % Compr + H-Matrix * Compr --> Compr
    elseif (Mh.typ==1) && (Ml.typ==0) && (Mr.typ==1)
        A  = alpha * (Ml * Mr.dat{1});
        B  = Mr.dat{2};
        Mh = hmxPlusAB(Mh,A,B);

    % Compr + H-Matrix * Full --> Full
    elseif (Mh.typ==1) && (Ml.typ==0) && (Mr.typ==2)
        error('hmxPlusMtimes : unvailable case')
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    % Compr + Compr * H-Matrix --> Compr
    elseif (Mh.typ==1) && (Ml.typ==1) && (Mr.typ==0)
        A  = alpha * Ml.dat{1}; 
        B  = Ml.dat{2} * Mr;
        Mh = hmxPlusAB(Mh,A,B);

    % Compr + Compr * Compr --> Compr
    elseif (Mh.typ==1) && (Ml.typ==1) && (Mr.typ==1) 
        A  = alpha * Ml.dat{1};
        B  = (Ml.dat{2} * Mr.dat{1}) * Mr.dat{2};
        Mh = hmxPlusAB(Mh,A,B);

    % Compr + Compr * Full --> Compr
    elseif (Mh.typ==1) && (Ml.typ==1) && (Mr.typ==2)
        A  = alpha * Ml.dat{1}; 
        B  = Ml.dat{2} * Mr.dat;
        Mh = hmxPlusAB(Mh,A,B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    % Compr + Full * H-Matrix --> Full
    elseif (Mh.typ==1) && (Ml.typ==2) && (Mr.typ==0)
        error('hmxPlusMtimes : unvailable case')
        
    % Compr + Full * Compr --> Compr
    elseif (Mh.typ==1) && (Ml.typ==2) && (Mr.typ==1)
        A  = alpha * (Ml.dat * Mr.dat{1});
        B  = Mr.dat{2};
        Mh = hmxPlusAB(Mh,A,B);

    % Compr + Full * Full --> Full
    elseif (Mh.typ==1) && (Ml.typ==2) && (Mr.typ==2)
        Mh.dat = full(Mh) + alpha * (Ml.dat * Mr.dat);
        Mh.typ = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % Full + --- --> Full
    elseif (Mh.typ==2) 
        Mh.dat = full(Mh) + alpha * (full(Ml) * full(Mr));
        Mh.typ = 2;

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
    % Unknown type    
    else
        error('hmxPlusMtimes : unvailable case')
    end
    
    
%%% Unavailable    
else
    error('hmxPlusMtimes : unvailable case')
end
end