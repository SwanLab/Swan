%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
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
%|    #    |   FILE       : nrtHmxCompressorSingular.m                    |
%|    #    |   VERSION    : 0.52                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 01.01.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : SVD-based compressors for singular matrix     |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Precision
tol = 1e-3;

% Read singular matrix
load('singular.mat');
size(M)
rank(M)

% Standard SVD
try
    [U,S,V] = svd(M);
catch
    disp('SVD did not converge')
end

% H-Matrix SVD 
[~,~,flag] = hmxSVD(M,tol);
flag 

% H-Matrix RSVD 
[~,~,flag] = hmxRSVD(M,tol,64);
flag 

% QR recompression
[A,B] = hmxQRSVD(M,eye(64),tol);




disp('~~> Michto gypsilab !')


