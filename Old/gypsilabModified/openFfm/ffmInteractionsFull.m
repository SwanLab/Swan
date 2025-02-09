function MV = ffmInteractionsFull(X,Xbox,Y,Ybox,V,Ibox,green,k)
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
%|    #    |   FILE       : ffmInteractionsFull.m                         |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : FUll product for non compressible leaves      |
%|  `---'  |                                                              |
%+========================================================================+

% Produit Matrice-Vecteur
MV = zeros(size(X,1),1,'like',V);

% Calcul plein pour chaque interaction
for i = 1:size(Ibox,1)
    % Boite en X 
    ix = Xbox.ind{Ibox(i,1)};
    nx = length(ix);
    
    % Boite en Y
    iy = Ybox.ind{Ibox(i,2)};
    ny = length(iy);
        
    % Noyau de green
    [idx,jdx] = ndgrid(ix',iy');
    Gxy       = ffmGreenKernel(X(idx(:),:),Y(jdx(:),:),green,k);
    Gxy       = reshape(Gxy,nx,ny);
    
    % Produit Matrice-Vecteur 
    MV(ix) = MV(ix) + reshape(Gxy,nx,ny) * V(iy);
end
end
