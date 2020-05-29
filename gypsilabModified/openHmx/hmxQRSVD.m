function [A,B] = hmxQRSVD(A,B,tol)
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
%|    #    |   FILE       : hmxQRSVD.m                                    |
%|    #    |   VERSION    : 0.52                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.01.2019                                    |
%| ( === ) |   SYNOPSIS   : QR factorization and SVD recompression for    |
%|  `---'  |                low-rank matrices                             |
%+========================================================================+

if ~isempty(A)
    % QR Factorisation A = QA * RA;
    [QA,RA] = qr(A,0);
    
    % QR Factorisation Bt = QB * RB
    [QB,RB] = qr(B.',0);
    
    % Singular value decomposition U*S*V = RA * RB.'
    try
        [U,S,V] = svd(RA * RB.','econ');
    catch
        return
    end

    % Compression rank
    s = diag(S);
    n = sum( s./s(1) >= tol );

    % Recompression A = QA * U * S
    A = QA * (U(:,1:n) * S(1:n,1:n));
    
    % Recompression B = V' * QB^t
    B = V(:,1:n)' * QB.';
end
end
