function Mh = hmxHorzcat(varargin)
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
%|    #    |   FILE       : hmxHorzcat.m                                  |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 15.07.2018                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix horizontal concatenation             |
%|  `---'  |                                                              |
%+========================================================================+

%%% Input analysis
Ml = varargin{1};
Mr = varargin{2};
if (nargin == 2)
    if isequal(Ml.pos{1},Mr.pos{1})
        Ix = (1:size(Ml,1))';        
        Iy = (1:size(Ml,2)+size(Mr,2))';
    else
        error('hmxHorzcat.m : unavailable case');
    end
else
    Ml = varargin{1};
    Mr = varargin{2};
    Ix = varargin{3};
    Iy = varargin{4};
end

%%% Preparation
% Particles
X = Ml.pos{1}(Ix,:);
Y = [Ml.pos{2} ; Mr.pos{2}];
Y = Y(Iy,:);

% Accuracy
tol = max(Ml.tol,Mr.tol);

% Initialiation
Mh = hmx(X,Y,tol);

% Admissibility
[isfar,Xdim,Ydim] = hmxFar(Mh);

% Y Indices
bool = (Iy<=size(Ml,2));
Iyl  = Iy(bool);
Iyr  = Iy(~bool)-size(Ml,2);

% Compression for far distances
if isfar
    % Compression
    [Al,Bl] = lowrank(Ml,Ix,Iyl);
    [Ar,Br] = lowrank(Mr,Ix,Iyr);
    
    % Concatenation
    A = [Al,Ar];
    B = zeros(size(A,2),length(Iy));
    B(1:size(Al,2),bool)      = Bl;
    B(size(Al,2)+1:end,~bool) = Br;
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);
    
    % Validation
    flag = (numel(A)+numel(B)) <= prod(size(Mh));
    
else
    flag = 0;
end


%%% Compression
if flag
    % Type
    Mh.typ = 1;
    
    % Low-rank
    Mh.dat = {A,B};


%%%% Full or sparse for smallest box (stopping criterion)
elseif (min(size(Mh)) < 100)
    % Type
    Mh.typ = 2;
    
    % Full
    Mh.dat = zeros(size(Mh));
    
    % Fill data
    Mh.dat(:,bool)  = full(Ml,Ix,Iyl);
    Mh.dat(:,~bool) = full(Mr,Ix,Iyr);


%%% H-Matrix (recursion)
else
    % Type
    Mh.typ = 0;
    
    % Subdivision for X
    [I1,I2] = hmxSubdivide(X,Xdim);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Y,Ydim);
    Mh.col  = {I1 , I2 , I1 , I2};
    
    % H-Matrix (recursion)
    for i = 1:4
        % Coordinates
        Ir = Mh.row{i};
        Ic = Mh.col{i};

        % Recursion
        Mh.chd{i} = hmxHorzcat(Ml,Mr,Ix(Ir),Iy(Ic)); 
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    
end
end






% function Mh = hmxHorzcat(Ml,Mr)
% %+========================================================================+
% %|                                                                        |
% %|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
% %|           openHmx is part of the GYPSILAB toolbox for Matlab           |
% %|                                                                        |
% %| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
% %| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
% %| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
% %| LICENCE   : This program is free software, distributed in the hope that|
% %| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
% %| redistribute and/or modify it under the terms of the GNU General Public|
% %| License, as published by the Free Software Foundation (version 3 or    |
% %| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
% %| is available, please contact us to activate a "pay for remove" option. |
% %| CONTACT   : matthieu.aussal@polytechnique.edu                          |
% %| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
% %|                                                                        |
% %| Please acknowledge the gypsilab toolbox in programs or publications in |
% %| which you use it.                                                      |
% %|________________________________________________________________________|
% %|   '&`   |                                                              |
% %|    #    |   FILE       : hmxHorzcat.m                                  |
% %|    #    |   VERSION    : 0.40                                          |
% %|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
% %|  ( # )  |   CREATION   : 14.03.2017                                    |
% %|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
% %| ( === ) |   SYNOPSIS   : Concatenate H-Matrix along 2nd dimension      |
% %|  `---'  |                with the rule                                 |
% %|         |                H-Matrix > Full > Compr                       |
% %+========================================================================+
% 
% % Check dimensions
% if ~(isa(Ml,'hmx') && isa(Mr,'hmx'))
%     error('hmxHorzcat.m : unavailable case')
% elseif (size(Ml,1)~=size(Mr,1))
%     error('hmxHorzcat.m : Dimensions of matrices being concatenated are not consistent')
% elseif norm(Ml.pos{1}-Mr.pos{1},'inf') > 1e-12
%     error('hmxHorzcat.m : unavailable case')
% end
%     
% % Initialization
% Mh = hmx(Ml.pos{1},[Ml.pos{2};Mr.pos{2}],Ml.tol);
% 
% 
% %%% [H-Matrix & H-Matrix] -> H-Matrix
% if (Ml.typ == 0) && (Mr.typ == 0)
%     % Particles distance average
%     Yl1 = mean(Ml.chd{1}.pos{2},1);
%     Yl2 = mean(Ml.chd{2}.pos{2},1);
%     Yr1 = mean(Mr.chd{1}.pos{2},1);
%     Yr2 = mean(Mr.chd{2}.pos{2},1);
%     
%     % Concatenation pairs
%     D     = [norm(Yl1-Yr1) norm(Yl1-Yr2) norm(Yl2-Yr1) norm(Yl2-Yr2)];
%     [r,d] = min(D);
%     if (d == 1) || (d == 4)
%         I = [1 1 ; 2 2 ; 3 3 ; 4 4];
%     else
%         I = [1 2 ; 2 1 ; 3 4 ; 4 3];
%     end
%     
%     % Construction (recursion)
%     Mh.typ = 0;
%     for i = 1:4
%         Mh.chd{i} = hmxHorzcat(Ml.chd{I(i,1)},Mr.chd{I(i,2)});
%         Mh.row{i} = Ml.row{i};
%         Mh.col{i} = [Ml.col{I(i,1)};size(Ml,2)+Mr.col{I(i,2)}];
%     end    
%     
% %     % Construction (recursion)
% %     Mh.typ = 0;
% %     for i = 1:4
% %         Mh.chd{i} = hmxHorzcat(Ml.chd{i},Mr.chd{i});
% %         Mh.row{i} = Ml.row{i};
% %         Mh.col{i} = [Ml.col{i};size(Ml,2)+Mr.col{i}];
% %     end    
%     
%     
% %     % Split for too big full leaf 
% %     for i = 1:4
% %         if (Mh.chd{i}.typ == 2) && (min(size(Mh.chd{i})) > 100)
% %             tmp       = Mh.chd{i}
% %             Mh.chd{i} = hmxBuilder(tmp.pos{1},tmp.pos{2},tmp.dat,tmp.tol);
% %         end
% %     end
% %     
% %     % Fusion
% %     Mh = hmxFusion(Mh);    
%     
% % [H-Matrix & Compr] -> H-Matrix           
% elseif (Ml.typ == 0) && (Mr.typ == 1)
%     Mr = hmxSplit(Ml.row,Mr);
%     Mh = hmxHorzcat(Ml,Mr);
%  
% % [H-Matrix & Full] -> H-Matrix 
% elseif (Ml.typ == 0) && (Mr.typ == 2)
%     Mr = hmxSplit(Ml.row,Mr);
%     Mh = hmxHorzcat(Ml,Mr);
% 
% 
% %%% [Compr & H-matrix] -> H-Matrix
% elseif (Ml.typ == 1) && (Mr.typ == 0)
%     Ml = hmxSplit(Mr.row,Ml);
%     Mh = hmxHorzcat(Ml,Mr);    
%     
% % [Compr & Compr] -> Compr
% elseif (Ml.typ == 1) && (Mr.typ == 1)
%     % Concatenation
%     nl = size(Ml.dat{1},2);
%     nr = size(Mr.dat{1},2);
%     A  = [Ml.dat{1} , Mr.dat{1}];
%     B  = [Ml.dat{2} , zeros(nl,size(Mr,2)) ;
%         zeros(nr,size(Ml,2)) , Mr.dat{2}] ;
%     
%     % Recompression
%     [A,B] = hmxQRSVD(A,B,Ml.tol);
%     
%     % Update
%     Mh.typ = 1;
%     Mh.dat = {A,B};
%     
%  
% % [Compr & Full] -> Full
% elseif (Ml.typ == 1) && (Mr.typ == 2)
%     Mh.typ = 2;
%     Mh.dat = [Ml.dat{1}*Ml.dat{2} , Mr.dat];    
% 
%     
% %%% [Full & H-matrix] -> H-Matrix
% elseif (Ml.typ == 2) && (Mr.typ == 0)
%     Ml = hmxSplit(Mr.row,Ml);
%     Mh = hmxHorzcat(Ml,Mr); 
%         
% % [Full & Compr] -> Full                            
% elseif (Ml.typ == 2) && (Mr.typ == 1)
%     Mh.typ = 2;
%     Mh.dat = [Ml.dat , Mr.dat{1}*Mr.dat{2}]; 
%         
% % [Full & Full] -> Full
% elseif (Ml.typ == 2) && (Mr.typ == 2)
%     Mh.typ = 2;
%     Mh.dat = [Ml.dat,Mr.dat];       
%     
% else
%     error('hmxHorzcat.m : unavailable case')
% end
% end
% 
% 
% function Mh = hmxSplit(row,Mh)
% % Subdivision for X
% Mh.row = row;
% 
% % Subdivision for Y
% [~,~,Ydim] = hmxFar(Mh);
% [I1,I2]    = hmxSubdivide(Mh.pos{2},Ydim);
% Mh.col     = {I1 , I2 , I1 , I2};
% 
% % H-Matrix (recursion)
% for i = 1:4
%     % Initialization
%     Mh.chd{i} = hmx(Mh.pos{1}(Mh.row{i},:),Mh.pos{2}(Mh.col{i},:),Mh.tol);
%     
%     % H-Matrix
%     if (Mh.typ == 0)
%         error('hmxSplit.m : unavailable case')
% 
%     % Compressed leaf
%     elseif (Mh.typ == 1)
%         % Subdivision
%         A = Mh.dat{1}(Mh.row{i},:);
%         B = Mh.dat{2}(:,Mh.col{i});
%         
%         % Recompression
%         [A,B] = hmxQRSVD(A,B,Mh.tol);
%         
%         % Update
%         Mh.chd{i}.typ = 1;
%         Mh.chd{i}.dat = {A,B};
%         Mh.chd{i}.tol = Mh.tol;
%         
%     % Full leaf
%     elseif (Mh.typ == 2)
%         Mh.chd{i}.typ = 2;
%         Mh.chd{i}.dat = Mh.dat(Mh.row{i},Mh.col{i});
%         Mh.chd{i}.tol = Mh.tol;
%     end
% end
% 
% % Erase data
% Mh.typ = 0;
% Mh.dat = [];
% end
