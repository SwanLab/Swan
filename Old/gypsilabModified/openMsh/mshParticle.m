function [mesh,elt2prt] = mshParticle(mesh)
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
%|    #    |   FILE       : mshParticle.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Submesh to particles mesh                     |
%|  `---'  |                                                              |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshParticle.m : unavailable case')
       
% Edge, triangular, tetrahedral mesh
elseif (size(mesh.elt,2) <= 4)
    % All particles
    prt2vtx = (1:size(mesh.vtx,1))';
    elt2prt = mesh.elt;
    
    % Colours
    col             = zeros(size(prt2vtx,1),1);
    tmp             = mesh.col * ones(1,size(mesh.elt,2));
    col(elt2prt(:)) = tmp(:);
    
    % Identify boundary
    bnd      = ( accumarray(elt2prt(:),1).*col ~= accumarray(elt2prt(:),tmp(:)) );
    col(bnd) = -1;
    
% Unknown type
else
    error('mshParticle.m : unavailable case')
end

% Mesh format
mesh = msh(mesh.vtx,prt2vtx,col);
end
