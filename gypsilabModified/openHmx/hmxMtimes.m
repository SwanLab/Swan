function Mh = hmxMtimes(Ml,Mr)
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
%|    #    |   FILE       : hmxMtimes.m                                   |
%|    #    |   VERSION    : 0.51                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Product of H-Matrix with the rule             |
%|  `---'  |                Compr > Full > H-Matrix                       |
%+========================================================================+

% Check dimensions
if (size(Ml,2) ~= size(Mr,1))
    error('hmxMtimes.m : matrix dimensions must agree.')
end

%%% H-Matrix * H-Matrix --> H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    % Initialisation
    Mh = hmx(Ml.pos{1},Mr.pos{2},Ml.tol);
    
    % H-Matrix * H-Matrix --> H-Matrix   (recursion)
    if (Ml.typ==0) && (Mr.typ==0)
        % Type
        Mh.typ = 0;
        
        % Bloc product indices
        I = [1 2; 1 2; 3 4; 3 4];
        J = [1 3; 2 4; 1 3; 2 4];
        
        % Bloc product
        for i = 1:4
            Mh.chd{i} = hmxMtimes(Ml.chd{I(i,1)},Mr.chd{J(i,1)}) + ...
                hmxMtimes(Ml.chd{I(i,2)},Mr.chd{J(i,2)});
            Mh.row{i} = Ml.row{I(i,1)};
            Mh.col{i} = Mr.col{J(i,2)};
        end
        
        % Fusion
        Mh = hmxFusion(Mh);
        
    % H-Matrix * Compr --> Compr
    elseif (Ml.typ==0) && (Mr.typ==1)
        Mh.typ = 1;
        Mh.dat = {hmxMtimes(Ml,Mr.dat{1}) , Mr.dat{2}};
        
    % H-Matrix * Full --> Full
    elseif (Ml.typ==0) && (Mr.typ==2)
        Mh.typ = 2;
        Mh.dat = hmxMtimes(Ml,Mr.dat);
       
        
    % Compr * H-Matrix --> Compr
    elseif (Ml.typ==1) && (Mr.typ==0)
        Mh.typ = 1;
        Mh.dat = {Ml.dat{1} , hmxMtimes(Ml.dat{2},Mr)};
            
    % Compr * Compr --> Compr
    elseif (Ml.typ==1) && (Mr.typ==1) 
        Mh.typ = 1;
        Mh.dat = {Ml.dat{1} , (Ml.dat{2} * Mr.dat{1}) * Mr.dat{2}};
        
    % Compr * Full --> Compr
    elseif (Ml.typ==1) && (Mr.typ==2)
        Mh.typ = 1;        
        Mh.dat = {Ml.dat{1} , Ml.dat{2} * Mr.dat};
        
        
    % Full * H-Matrix --> Full
    elseif (Ml.typ==2) && (Mr.typ==0)
        Mh.typ = 2;
        Mh.dat = hmxMtimes(Ml.dat,Mr);
        
    % Full * Compr --> Compr
    elseif (Ml.typ==2) && (Mr.typ==1)
        Mh.typ = 1;
        Mh.dat = {Ml.dat*Mr.dat{1} , Mr.dat{2}};
        
    % Full * Full --> Full or Compr
    elseif (Ml.typ==2) && (Mr.typ==2)
        if (min(size(Mh)) < 100)
             Mh.typ = 2;
             Mh.dat = Ml.dat * Mr.dat;
        else
             Mh.typ = 1;
             [A,B]  = hmxQRSVD(Ml.dat,Mr.dat,Mh.tol);
             Mh.dat = {A,B};
        end
  
        
    % Unknown type    
    else
        error('hmxMtimes : unvailable case')
    end
    
    
%%% H-Matrix * Matrix --> Full
elseif isa(Ml,'hmx')
    if (size(Mr,2) > 0)
        % H-Matrix (recursion)
        if (Ml.typ == 0)
            % Initializaton
            Mh = zeros(size(Ml,1),size(Mr,2),class(Ml.row{1}));
            
            % Recursion
            for i = 1:4
                Mh(Ml.row{i},:) = Mh(Ml.row{i},:) + hmxMtimes(Ml.chd{i},Mr(Ml.col{i},:));
            end
            
        % Compressed leaf
        elseif (Ml.typ == 1)
            Mh = Ml.dat{1} * (Ml.dat{2} * Mr);
            
        % Full leaf
        elseif (Ml.typ == 2)
            Mh = Ml.dat * Mr;
            
        % Unknown type
        else
            error('hmxMtimes.m : unavailable case')
        end
    else
        Mh = zeros(size(Ml,1),0);
    end
    
    
%%% Matrix * H-Matrix --> Full
elseif isa(Mr,'hmx')
    if (size(Ml,1) > 0)
        % H-Matrix (recursion)
        if (Mr.typ == 0)
            % Initializaton
            Mh = zeros(size(Ml,1),size(Mr,2),class(Mr.row{1}));
            
            % Recursion
            for i = 1:4
                Mh(:,Mr.col{i}) = Mh(:,Mr.col{i}) + hmxMtimes(Ml(:,Mr.row{i}),Mr.chd{i});
            end
            
        % Compressed leaf
        elseif (Mr.typ == 1)
            Mh = (Ml * Mr.dat{1}) * Mr.dat{2};
            
        % Full leaf
        elseif (Mr.typ == 2)
            Mh = Ml * Mr.dat;
            
        % Unknown type
        else
            error('hmxMtimes.m : unavailable case')
        end
    else
        Mh = zeros(0,size(Mr,2));
    end
end
end
