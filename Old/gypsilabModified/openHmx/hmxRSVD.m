function [A,B,flag] = hmxRSVD(varargin)
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
%|    #    |   FILE       : hmxRSVD.m                                     |
%|    #    |   VERSION    : 0.52                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.01.2019                                    |
%| ( === ) |   SYNOPSIS   : RSVD from Halko et al. 'finding structure with|
%|  `---'  |                randomness'. Adapted from Antoine Liutkus.    |
%+========================================================================+

% Input analysis
if (numel(varargin{2}) == 1)
    % Input
    M   = varargin{1};
    tol = varargin{2};
    rk  = varargin{3};
    
    % Dimensions
    [m,n] = size(M);
    
    % Matrix vectot product
    MV = @(V) M * V;
    VM = @(V) V * M;
    
    % Type
    typ = class(M);
    
else
    % Input
    A   = varargin{1};
    B   = varargin{2};
    tol = varargin{3};
    
    % Dimensions
    [m,rk] = size(A);
    n      = size(B,2);
    
    % Type
    typ = class(A);
    
    % Matrix vectot product
    MV    = @(V) A * (B * V);
    VM    = @(V) (V * A) * B;
end

% Randomized
p = min(2*rk,n);
X = randn(n,p,typ);
Y = MV(X);
try
    W1 = orth(Y);
catch
    A    = [];
    B    = [];
    flag = 0;
    return
end
B = VM(W1');

% Truncated SVD
try
    [W2,S,V] = svd(B,'econ');
catch
    A    = [];
    B    = [];
    flag = 0;
    return
end

% Product 
U = W1*W2;

% Rank with fixed accuracy 
if (numel(S) ~= 0)
    I = find(abs(diag(S)/S(1)) >= tol);
else
    I = [];
end
    
% No values
if isempty(I)
   A    = zeros(m,0);
   B    = zeros(0,n);
   flag = 1;

% Accuracy not reached
elseif (length(I) >= rk)
    A    = [];
    B    = [];
    flag = 0;
    
% Low-rank representation
else
    A    = U(:,I);
    B    = S(I,I) * V(:,I)';
    flag = 1;
end
end
