classdef ffm
%+========================================================================+
%|                                                                        |
%|         OPENFFM - LIBRARY FOR FAST AND FREE MEMORY CONVOLUTION         |
%|           openFfm is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : ffm.m                                         |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Fast & Furious Method object                  |
%|  `---'  |                                                              |
%+========================================================================+

properties
    dim = [];          % MATRIX DIMENSION
    fcn = @(V) V;      % MATRIX-VECTOR PRODUCT FUNCTION
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Mv = ffm(varargin)
        % Empty initialization
        if (nargin == 0)
            
        % Particle builder
        elseif (nargin == 5)
            X      = varargin{1};
            Y      = varargin{2};
            green  = varargin{3};
            k      = varargin{4};
            tol    = varargin{5};
            Mv.dim = [size(X,1) size(Y,1)];
            Mv.fcn = @(V) ffmProduct(X,Y,V,green,k,tol);

        % Finite element builder
        elseif (nargin == 7)
            Mx     = varargin{1};
            X      = varargin{2};
            green  = varargin{3};
            k      = varargin{4};            
            Y      = varargin{5};
            My     = varargin{6};
            tol    = varargin{7};
            Mv.dim = [size(Mx,1) size(My,2)];
            Mv.fcn = @(V) Mx*ffmProduct(X,Y,My*V,green,k,tol);

        else
            error('ffm.m : undefined constructor case')
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIZE
    function s = size(varargin)
        s = varargin{1}.dim;
        if (nargin == 2)
            s = s(varargin{2});
        end
    end
    
    % LENGTH
    function l = length(Mv)
        l = max(size(Mv));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGEBRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDITION
    function Ml = plus(Ml,Mr)
        if isa(Ml,'ffm') && isa(Mr,'ffm')
            Ml.fcn = @(V) Ml.fcn(V) + Mr.fcn(V);
        elseif isa(Ml,'ffm') && isnumeric(Mr)
            Ml.fcn = @(V) Ml.fcn(V) + Mr*V;
        elseif isnumeric(Ml) && isa(Mr,'ffm')
            Mr.fcn = @(V) Ml*V + Mr.fcn(V);
            Ml     = Mr;
        else
            error('ffm.m : unavailable case')
        end
    end
    
    % SUBSTRACTION
    function Ml = minus(Ml,Mr)
        Ml = Ml + (-Mr);
    end    
    
    % SCALAR PRODUCT
    function Ml = times(Ml,Mr)
        if isa(Ml,'ffm') && isa(Mr,'ffm')
            error('ffm.m : unavailable case')
        elseif isa(Ml,'ffm') && isnumeric(Mr) && (numel(Mr)==1)
            Ml.fcn = @(V) Ml.fcn(Mr.*V);
        elseif isnumeric(Ml) && (numel(Ml)==1) && isa(Mr,'ffm')
            Mr.fcn = @(V) Ml.*Mr.fcn(V);
            Ml     = Mr;
        else
            error('ffm.m : unavailable case')
        end
    end
    
    % UMINUS 
    function Mv = uminus(Mv)
        Mv = times(-1,Mv);
    end
    
    % MATRIX PRODUCT
    function Ml = mtimes(Ml,Mr)
        if isa(Ml,'ffm') && isa(Mr,'ffm')
            Ml.dim = [Ml.dim(1) Mr.dim(2)];
            Ml.fcn = @(V) Ml.fcn(Mr.fcn(V));
            
        elseif isa(Ml,'ffm') && issparse(Mr)
            Ml.dim = [Ml.dim(1) size(Mr,2)];
            Ml.fcn = @(V) Ml.fcn(Mr*V);

        elseif isa(Ml,'ffm') && isnumeric(Mr) && (numel(Mr)==1)
            Ml = times(Ml,Mr);
            
        elseif isa(Ml,'ffm') && isnumeric(Mr)
            MV = zeros(Ml.dim(1),size(Mr,2));
            for i = 1:size(MV,2)
                MV(:,i) = Ml.fcn(Mr(:,i));
            end
            Ml = MV;
            
        elseif isnumeric(Ml) && (numel(Ml)==1) && isa(Mr,'ffm')
            Ml = times(Ml,Mr);
            
        elseif isnumeric(Ml) && isa(Mr,'ffm') 
            Mr.dim = [size(Ml,1) Mr.dim(2)];
            Mr.fcn = @(V) Ml*Mr.fcn(V);
            Ml     = Mr;            
        else
            error('ffm.m : unavailable case')
        end
    end
    
end
end
