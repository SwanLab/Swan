function [A,B] = hmxLowrank(Mh)
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
%|    #    |   FILE       : hmxLowrank.m                                  |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix to low-rank approximation    |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix (recursion)
if (Mh.typ == 0)
    % Initialisation
    A = zeros(size(Mh,1),0,class(Mh.row{1}));
    B = zeros(0,size(Mh,2),class(Mh.row{1}));

    % Recursion
    for i = 1:4
        % Children evaluation
        [Ai,Bi] = hmxLowrank(Mh.chd{i});

        % Low-rank addition for A
        tmp              = zeros(size(Mh,1),size(Ai,2),class(Mh.row{1}));
        tmp(Mh.row{i},:) = Ai;
        A                = [A , tmp];

        % Low-rank addition for B
        tmp              = zeros(size(Bi,1),size(Mh,2),class(Mh.row{1}));
        tmp(:,Mh.col{i}) = Bi;
        B                = [B ; tmp];
    end

    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);

%%% Compressed leaf
elseif (Mh.typ == 1)
    A = Mh.dat{1};
    B = Mh.dat{2};

%%% Full leaf
elseif (Mh.typ == 2)
    % Low-Rank conversion
    A = full(Mh.dat);
    B = eye(size(Mh,2));
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);

%%% Others    
else
    error('hmxLowrank.m : unavailable case')
end
end
