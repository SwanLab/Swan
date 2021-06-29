function [mesh1,mesh2] = mshSplit(mesh,X0,U)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : mshSplit.m                                    |
%|    #    |   VERSION    : 0.51                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 21.06.2019                                    |
%| ( === ) |   SYNOPSIS   : Split mesh with planar cut                    |
%|  `---'  |                (to be clearly improved)                      |
%+========================================================================+

% Normalize direction
N = U./norm(U);

% Split mesh
lambda = (mesh.ctr-ones(length(mesh),1)*X0) * N';
ind    = (lambda>0);
mesh1  = mesh.sub(ind);
mesh2  = mesh.sub(~ind);

% Extract free boundary
bound  = mesh1.bnd;
mu     = (bound.vtx-ones(size(bound.vtx,1),1)*X0) * N';

% Move upper boundary
[~,IA] = intersect(msh(mesh1.vtx),msh(bound.vtx));
mesh1.vtx(IA,:) = mesh1.vtx(IA,:) - mu*N;

% Move lower boundary
[~,IA] = intersect(msh(mesh2.vtx),msh(bound.vtx));
mesh2.vtx(IA,:) = mesh2.vtx(IA,:) - mu*N;
end
