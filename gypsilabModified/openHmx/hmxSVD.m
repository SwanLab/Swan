function [A,B,flag] = hmxSVD(M,tol)
%%+========================================================================+
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
%|    #    |   FILE       : hmxSVD.m                                      |
%|    #    |   VERSION    : 0.52                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.01.2019                                    |
%| ( === ) |   SYNOPSIS   : SVD compression for full matrix               |
%|  `---'  |                                                              |
%+========================================================================+

% Empty matrix
if isempty(M)
    A = [];
    B = [];

else
    % Singular value decomposition
    try
        [U,S,V] = svd(M,'econ');
    catch
        A    = [];
        B    = [];
        flag = 0;
        return
    end
    
    % Rank estimation at machine precision
    s   = diag(S);
    acc = max(size(M)) * eps(max(s));
    rk  = sum(s > acc);
    
    % No rank default
    if (rk == length(s))
        A    = [];
        B    = [];
        flag = 0;
        
    % Rank default
    else
        % Compression rank
        n = sum( s./s(1) >= tol );
        
        % Precision not reached
        if (n*sum(size(M)) > numel(M))
            A    = [];
            B    = [];
            flag = 0;
            
        % Compression
        else
            A    = U(:,1:n);
            B    = S(1:n,1:n)*V(:,1:n)';
            flag = 1;
        end
    end
end
