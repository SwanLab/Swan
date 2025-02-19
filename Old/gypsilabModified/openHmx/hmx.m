classdef hmx
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
%|    #    |   FILE       : hmx.m                                         |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix class definition and functions       |
%|  `---'  |                                                              |
%+========================================================================+

properties
    typ = [];            % LEAF TYPE (0=H-MATRIX ; 1=COMPRESSED ; 2=FULL/SPARSE) 
    pos = [];            % COORDINATES POSITIONS (X,Y)
    row = cell(1,4);     % CHILDREN ROWS INDICES (M11,M12,M21,M22)
    col = cell(1,4);     % CHILDREN COLUMNS INDICES
    chd = cell(1,4);     % CHILDREN 
    dat = [];            % LEAF DATA
    tol = [];            % COMPRESSORS ACCURACY
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONSTRUCTOR
    function Mh = hmx(varargin)
        % Empty object
        if (nargin == 0)
            
        % Initialization other H-matrix
        elseif (nargin == 2)
            Mh = hmxCopy(varargin{1},varargin{2});
            
        % Initialization with dimension and accuracy    
        elseif (nargin == 3)
            Mh     = hmx();
            Mh.pos = {varargin{1},varargin{2}};
            Mh.tol = varargin{3};
       
        % Particles builder with partial and total pivoting   
        elseif (nargin == 4)
            X     = varargin{1};
            Y     = varargin{2};
            green = varargin{3};
            acc   = varargin{4};
            Mh    = hmxBuilder(X,Y,green,acc);
            
        % Compressed builder    
        elseif (nargin == 5)
            X      = varargin{1};
            Y      = varargin{2};
            A      = varargin{3};
            B      = varargin{4};
            acc    = varargin{5};
            Mh     = hmx(X,Y,acc);
            Mh.typ = 1;
            Mh.dat = {A,B};            
            
       % Finite element builder    
        elseif (nargin == 8)
            Xunk  = varargin{1};
            Yunk  = varargin{2};
            Mx    = varargin{3};
            X     = varargin{4};
            green = varargin{5};
            Y     = varargin{6};
            My    = varargin{7};
            acc   = varargin{8};
            Mh    = hmxBuilderFem(Xunk,Yunk,Mx,X,green,Y,My,acc);
            
        else
            error('hmx.m : undefined constructor case')
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE
    function s = size(varargin)
        Mh = varargin{1};
        s  = [size(Mh.pos{1},1),size(Mh.pos{2},1)];
        if (nargin == 2)
            s = s(varargin{2});
        end
    end
    
    % LENGTH
    function l = length(Mh)
        l = max(size(Mh));
    end

    % ISLOWER
    function b = islower(Mh)
        if (Mh.typ == 0)
            b = (Mh.chd{2}.typ == 1) && isempty(Mh.chd{2}.dat{1});
        elseif (Mh.typ == 2)
            b = 1;
        else
            error('hmx.m : unavailable case.')
        end
    end
        
    % ISUPPER
    function b = isupper(Mh)
        if (Mh.typ == 0)
            b = (Mh.chd{3}.typ == 1) && isempty(Mh.chd{3}.dat{1});
        elseif (Mh.typ == 2)
            b = 1;
        else
            error('hmx.m : unavailable case.')
        end
    end    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STRUCTURE VISUALISATION
    function spy(Mh)
        hmxSpy(Mh);
    end
    
    % PLOT POSITIONS
    function plot3(Mh)
        plot3(Mh.pos{1}(:,1),Mh.pos{1}(:,2),Mh.pos{1}(:,3),'bo',...
            Mh.pos{2}(:,1),Mh.pos{2}(:,2),Mh.pos{2}(:,3),'*r');
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FULL
    function M = full(varargin)
        if (nargin == 1)
            M = hmxFull(varargin{1});
        else
            M = hmxFullSub(varargin{1},varargin{2},varargin{3});
        end
    end
    
    % SPARSE
    function M = sparse(varargin)
        if (nargin == 1)
            M = hmxSparse(varargin{1});
        elseif (nargin == 2)
            M = hmxSparsify(varargin{1},varargin{2});
        else
            M = hmxSparseSub(varargin{1},varargin{2},varargin{3});
        end
    end
    
    % LOW-RANK
    function [A,B] = lowrank(varargin)
        if (nargin == 1)
            [A,B] = hmxLowrank(varargin{1});
        else
            [A,B] = hmxLowrankSub(varargin{1},varargin{2},varargin{3});
        end
    end
    
    % DIAGONAL
    function D = diag(Mh)
        D = hmxDiag(Mh,(1:size(Mh,1))',(1:size(Mh,2))');
        D = full(diag(D));
    end
    
    % ZEROS
    function Mh = zeros(Mh)
        Mh.typ = 1;
        Mh.row = cell(1,4);
        Mh.col = cell(1,4);
        Mh.chd = cell(1,4);
        Mh.dat = {zeros(size(Mh,1),0),zeros(0,size(Mh,2))};
    end
    
    % ONES
    function Mh = ones(Mh)
        Mh.typ = 1;
        Mh.row = cell(1,4);
        Mh.col = cell(1,4);
        Mh.chd = cell(1,4);
        Mh.dat = {ones(size(Mh,1),1),ones(1,size(Mh,2))};
    end

    % SINGLE
    function Mh = single(Mh)
        Mh = hmxSingle(Mh);
    end

    % DOUBLE
    function Mh = double(Mh)
        Mh = hmxDouble(Mh);
    end

    % TRANSPOSITION
    function Mh = transpose(Mh)
        Mh = hmxTranspose(Mh);
    end
    
    % TRANSPOSITION CONJUGATE
    function Mh = ctranspose(Mh)
        Mh = hmxCtranspose(Mh);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGEBRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCALAR PRODUCT
    function Mh = times(Ml,Mr)
        Mh = hmxTimes(Ml,Mr);
    end
    
    % UMINUS 
    function Mh = uminus(Mh)
        Mh = (-1).*Mh;
    end
    
    % ADDITION
    function Mh = plus(Ml,Mr)
        Mh = hmxPlus(Ml,Mr);
    end
    
    % SUBSTRACTION
    function Mh = minus(Ml,Mr)
       Mh = Ml + (-Mr); 
    end
    
    % MATRIX PRODUCT
    function Mh = mtimes(Ml,Mr)
        if isnumeric(Ml) && (numel(Ml) == 1)
            Mh = Ml .* Mr;
        elseif isnumeric(Mr) && (numel(Mr) == 1)
            Mh = Ml .* Mr;
        else
            Mh = hmxMtimes(Ml,Mr);
        end
    end
    
    % ADDITION WITH MATRIX PRODUCT
    function Mh = plusmtimes(Mh,alpha,Ml,Mr)
        Mh = hmxPlusMtimes(Mh,alpha,Ml,Mr);
    end
    
    % INVERSION
    function Mh = inv(Mh)
        Mh = hmxInv(Mh);
    end
    
    % CHOLESKY FACTORISATION
    function Mh = chol(Mh)
        Mh = hmxChol(Mh);
    end
    
    % LDLt FACTORISATION
    function [Mh,Dh] = ldl(Mh)
        [Mh,Dh] = hmxLdl(Mh);
    end
    
    % LU FACTORISATION
    function [Lh,Uh] = lu(Mh)
        [Lh,Uh] = hmxLU(Mh);
    end
    
    % MLDIVIDE
    function B = mldivide(Mh,B)
       B = hmxMldivide(Mh,B); 
    end
        
    % MRDIVIDE
    function B = mrdivide(B,Mh)
        B = (Mh.'\B.').';
    end  
    
    % VERTICAL CONCATENATION
    function Mh = vertcat(varargin)
        Mh = varargin{1}.'; 
        for i = 2:nargin
            Mh = hmxHorzcat(Mh,varargin{i}.');
        end
        Mh = Mh.';
    end
    
    % HORIZONTAL CONCATENATION
    function Mh = horzcat(varargin)
        Mh = varargin{1}; 
        for i = 2:nargin
            Mh = hmxHorzcat(Mh,varargin{i});
        end
    end
end
end
