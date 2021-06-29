classdef bmm
%+========================================================================+
%|                                                                        |
%|               OPENBMM - LIBRARY FOR BLOCK MATRIX ALGEBRA               |
%|           openBmm is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : bmm.m                                         |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Block matrix Method object                    |
%|  `---'  |                                                              |
%+========================================================================+

properties
    row = {};      % ROWS INDICES
    col = {};      % COLUMNS INDICES
    blk = {};      % BLOCK DATA
end

methods  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Mb = bmm(varargin)
        % Empty object
        if (nargin == 0)
                        
        % Initialization with deep copy, single block and multi-block 
        elseif (nargin == 1)
            Mb = bmmBuilder(varargin{1});

        % Initialization with indices and block
        elseif (nargin == 3) || (nargin==4)
            Mb     = bmm();
            Mb.row = varargin{1};
            Mb.col = varargin{2};
            if iscell(varargin{3})
                Mb.blk = varargin{3};
            else
                Mb.blk = cell(length(Mb.row),length(Mb.col));
                if isa(varargin{3},'Composite')
                    for i = 1:length(Mb.row)
                        if iscell(varargin{3}{i}) && (nargin==4)
                            tmp = varargin{3}{i};
                            for j = 1:length(tmp)
                                Mb.blk{i,j} = hmx(tmp{j},varargin{4});
                            end
                        else
                            Mb.blk(i,:) = varargin{3}{i};
                        end
                    end
                else
                    for i = 1:length(Mb.row)
                        for j = 1:length(Mb.col)
                            Mb.blk{i,j} = varargin{3}(Mb.row{i},Mb.col{j});
                        end
                    end
                end
            end
            
        else
            error('bmm.m : undefined constructor case')
        end
    end
    
    % ZEROS
    function Mb = zeros(Mb)
        for n = 1:numel(Mb.blk)
            if isa(Mb.blk{n},'hmx')
                Mb.blk{n} = zeros(Mb.blk{n});
            else
                Mb.blk{n} = sparse(size(Mb.blk{n},1),size(Mb.blk{n},2));
            end
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE
    function s = size(varargin)
        I = cell2mat(varargin{1}.row);
        J = cell2mat(varargin{1}.col);
        s = [length(I) length(J)];
        if (nargin == 2)
            s = s(varargin{2});
        end
    end
    
    % LENGTH
    function l = length(Mb)
        l = max(size(Mb));
    end

    % CLASS
    function cls = class(Mb)
        cls = cell(size(Mb.blk));
        for n = 1:numel(Mb.blk)
            cls{n} = class(Mb.blk{n});
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STRUCTURE VISUALISATION
    function spy(Mb)
        bmmSpy(Mb);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FULL CONVERSION
    function M = full(Mb)
        M = zeros(size(Mb));
        for i = 1:length(Mb.row)
            for j = 1:length(Mb.col)
                M(Mb.row{i},Mb.col{j}) = full(Mb.blk{i,j});
            end
        end
    end
    
    % SPARSE CONVERSION
    function M = sparse(Mb)
        M = sparse(size(Mb));
        for i = 1:length(Mb.row)
            for j = 1:length(Mb.col)
                M(Mb.row{i},Mb.col{j}) = sparse(Mb.blk{i,j});
            end
        end
    end

    % TRANSPOSITION
    function Mb = transpose(Mb)
        tmp    = Mb.row;
        Mb.row = Mb.col;
        Mb.col = tmp;
        Mb.blk = Mb.blk';
        for n = 1:numel(Mb.blk)
            Mb.blk{n} = Mb.blk{n}.';
        end
    end
    
    % TRANSPOSITION CONJUGATE
    function Mb = ctranspose(Mb)
        tmp    = Mb.row;
        Mb.row = Mb.col;
        Mb.col = tmp;
        Mb.blk = Mb.blk';
        for n = 1:numel(Mb.blk)
            Mb.blk{n} = Mb.blk{n}';
        end
    end
    
    % VERTICAL CONCATENATION
    function Mb = vertcat(varargin)
        Mb = varargin{1}; 
        for i = 2:nargin
            Mb = cat(1,Mb,varargin{i});
        end
    end
    
    % HORIZONTAL CONCATENATION
    function Mb = horzcat(varargin)
        Mb = varargin{1}; 
        for i = 2:nargin
            Mb = cat(2,Mb,varargin{i});
        end
    end
    
    % CONCATENATION
    function Ml = cat(dim,Ml,Mr)
        Ml = bmmCat(dim,Ml,Mr);
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGEBRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCALAR PRODUCT
    function Ml = times(Ml,Mr)
        Ml = bmmTimes(Ml,Mr);
    end
    
    % UMINUS
    function Mb = uminus(Mb)
        Mb = times(Mb,-1);
    end
    
    % ADDITION
    function Ml = plus(Ml,Mr)
        Ml = bmmPlus(Ml,Mr);
    end
    
    % SUBSTRACTION
    function Ml = minus(Ml,Mr)
        Ml = Ml + (-Mr);
    end
    
    % MATRIX PRODUCT
    function M = mtimes(Ml,Mr)
        if isnumeric(Ml) && (numel(Ml) == 1)
            M = bmmTimes(Mr,Ml);
        elseif isnumeric(Mr) && (numel(Mr) == 1)
            M = bmmTimes(Ml,Mr);
        else
            M = bmmMtimes(Ml,Mr);
        end
    end
    
    % INVERSION
    function Mb = inv(Mb)
        Ib = bmm(Mb.row,Mb.col,speye(size(Mb)));
        for n = 1:numel(Mb.blk)
            if isa(Mb.blk{n},'hmx')
                Ib.blk{n} = hmx(Mb.blk{n}.pos{1},Mb.blk{n}.pos{2},Ib.blk{n},Mb.blk{n}.tol);
            elseif isa(Mb.blk,'fmm')
                error('bmm.m : unavailable case')
            end
        end
        Mb = Mb\Ib;
    end
    
    % LU FACTORISATION
    function [Lb,Ub] = lu(Mb)
        [Lb,Ub] = bmmLU(Mb);
    end
    
    % MLDIVIDE
    function B = mldivide(Mb,B)
        B = bmmMldivide(Mb,B);
    end
end
end
