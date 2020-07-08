function Bh = hmxMldivide(Mh,Bh)
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
%|    #    |   FILE       : hmxMldivide.m                                 |
%|    #    |   VERSION    : 0.51                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.01.2019                                    |
%| ( === ) |   SYNOPSIS   : Solve Mh x = Bh with the rule                 |
%|  `---'  |                Compr > Full > H-Matrix                       |
%+========================================================================+

% Check dimensions
if (size(Mh,1) ~= size(Bh,1))
    error('hmxMldivide.m : matrix dimensions must agree.')
end

% Check input data
if ~isa(Mh,'hmx')
    error('hmxMldivide.m : unavailable case.')
end

%%% Lower H-Matrix \ H-Matrix -> H-Matrix
if islower(Mh) && isa(Bh,'hmx')  
    % H-Matrix \ H-Matrix -> H-Matrix (recursion)
    if (Mh.typ == 0) && (Bh.typ == 0) 
        % X11 -> L11 \ B11
        Bh.chd{1} = hmxMldivide(Mh.chd{1},Bh.chd{1});
        
        % X12 -> L11 \ B12
        Bh.chd{2} = hmxMldivide(Mh.chd{1},Bh.chd{2});
        
        % X21 -> L22 \ (B21 - L21*X11)
        Bh.chd{3} = Bh.chd{3} - Mh.chd{3} * Bh.chd{1};
        Bh.chd{3} = hmxMldivide(Mh.chd{4},Bh.chd{3});
        
        % X22 -> L22 \ (B22 - L21 * X12)
        Bh.chd{4} = Bh.chd{4} - Mh.chd{3} * Bh.chd{2};
        Bh.chd{4} = hmxMldivide(Mh.chd{4},Bh.chd{4}); 
        
        % Fusion
        Bh = hmxFusion(Bh);
        
    % H-Matrix \ Compr -> Compr  
    elseif (Mh.typ == 0) && (Bh.typ == 1) 
        Bh.dat = { hmxMldivide(Mh,Bh.dat{1}) , Bh.dat{2} };
        
    % H-Matrix \ Full -> Full
    elseif (Mh.typ == 0) && (Bh.typ == 2) 
        Bh.typ = 2;
        Bh.dat = hmxMldivide(Mh,Bh.dat);
        

    % Compr \ --- -> ---
    elseif (Mh.typ == 1) 
        error('hmxMldivide : unavailable case')
                
        
    % Full \ H-Matrix -> Unavailable
    elseif (Mh.typ == 2) && (Bh.typ == 0)
        error('hmxMldivide : unavailable case')
        
    % Full \ Compr -> Compr
    elseif (Mh.typ == 2) && (Bh.typ == 1)
        Bh.dat = { Mh.dat\Bh.dat{1} , Bh.dat{2} };
        
    % Full \ Full -> Full
    elseif (Bh.typ == 2)
        Bh.dat = Mh.dat \ Bh.dat;
        
        
    else
        error('hmxMldivide : unavailable case')
    end


%%% Upper H-Matrix \ H-Matrix -> H-Matrix
elseif isupper(Mh) && isa(Bh,'hmx')  
    % H-Matrix \ H-Matrix -> H-Matrix (recursion)
    if (Mh.typ == 0) && (Bh.typ == 0) 
        % X22 -> U22 \ B22
        Bh.chd{4} = hmxMldivide(Mh.chd{4},Bh.chd{4});
        
        % X21 -> U22 \ B21
        Bh.chd{3} = hmxMldivide(Mh.chd{4},Bh.chd{3});
        
        % X12 -> U11 \ (B12 - U12*X22)
        Bh.chd{2} = Bh.chd{2} - Mh.chd{2} * Bh.chd{4};
        Bh.chd{2} = hmxMldivide(Mh.chd{1},Bh.chd{2});
        
        % X11 -> U11 \ (B11 - U12 * X21)
        Bh.chd{1} = Bh.chd{1} - Mh.chd{2} * Bh.chd{3};
        Bh.chd{1} = hmxMldivide(Mh.chd{1},Bh.chd{1}); 
        
        % Fusion
        Bh = hmxFusion(Bh);
        
    % H-Matrix \ Compr -> Compr  
    elseif (Mh.typ == 0) && (Bh.typ == 1) 
        Bh.dat = { hmxMldivide(Mh,Bh.dat{1}) , Bh.dat{2} };
        
    % H-Matrix \ Full -> Full
    elseif (Mh.typ == 0) && (Bh.typ == 2) 
        Bh.typ = 2;
        Bh.dat = hmxMldivide(Mh,Bh.dat);
        

    % Compr \ --- -> ---
    elseif (Mh.typ == 1) 
        error('hmxMldivide : unavailable case')
                
        
    % Full \ H-Matrix -> Unavailable
    elseif (Mh.typ == 2) && (Bh.typ == 0)
        error('hmxMldivide : unavailable case')
        
    % Full \ Compr -> Compr
    elseif (Mh.typ == 2) && (Bh.typ == 1)
        Bh.dat = { Mh.dat\Bh.dat{1} , Bh.dat{2} };
        
    % Full \ Full -> Full
    elseif (Bh.typ == 2)
        Bh.dat = Mh.dat \ Bh.dat;
        
        
    else
        error('hmxMldivide : unavailable case')
    end
    
    
%%% Lower H-Matrix \ Matrix -> Matrix   
elseif islower(Mh)
    if (size(Bh,2) > 0)
        % H-Matrix (recursion)
        if (Mh.typ == 0)
            % X1 -> L11 \ B1
            X1 = hmxMldivide(Mh.chd{1},Bh(Mh.row{1},:));
            
            % X2 -> L22 \ (B2 - L21*X1)
            X2 = Bh(Mh.row{3},:) - Mh.chd{3}*X1;
            X2 = hmxMldivide(Mh.chd{4},X2);
            
            % Bh = [X1 X2]
            if issparse(X1)
                Bh = sparse(size(Mh,1),size(Bh,2));
            else
                Bh = zeros(size(Mh,1),size(Bh,2),class(X1));
            end
            Bh(Mh.col{1},:) = X1;
            Bh(Mh.col{2},:) = X2;
            
        % Compressed leaf
        elseif (Mh.typ == 1)
            error('hmxMldivide : unavailable case')
            
        % Full leaf
        elseif (Mh.typ == 2)
            Bh = Mh.dat \ Bh;
            
        % Unknown type
        else
            error('hmxMldivide : unavailable case')
        end
    else
        Bh = zeros(size(Mh,1),0);
    end
    
    
%%% Upper H-Matrix \ Matrix -> Matrix   
elseif isupper(Mh)
    if (size(Bh,2) > 0)
        % H-Matrix (recursion)
        if (Mh.typ == 0)
            % X2 -> U22 \ B2
            X2 = hmxMldivide(Mh.chd{4},Bh(Mh.row{4},:));
            
            % X1 -> U11 \ (B1 - U12*X2)
            X1 = Bh(Mh.row{1},:) - Mh.chd{2}*X2;
            X1 = hmxMldivide(Mh.chd{1},X1);
            
            % Bh = [X1 ; X2]
            if issparse(X1)
                Bh = sparse(size(Mh,1),size(Bh,2));
            else
                Bh = zeros(size(Mh,1),size(Bh,2),class(X1));
            end
            Bh(Mh.col{1},:) = X1;
            Bh(Mh.col{2},:) = X2;
            
        % Compressed leaf
        elseif (Mh.typ == 1)
            error('hmxMldivide : unavailable case')
            
        % Full leaf
        elseif (Mh.typ == 2)
            Bh = Mh.dat \ Bh;
            
        % Unknown type
        else
            error('hmxMldivide : unavailable case')
        end
    else
        Bh = zeros(size(Mh,1),0);
    end
    
    
%%% H-Matrix \ (H-)Matrix -> (H-)Matrix
else
    [Lh,Uh] = lu(Mh);
    Bh      = Uh\(Lh\Bh);
end
