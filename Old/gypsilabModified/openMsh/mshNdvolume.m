function V = mshNdvolume(mesh)
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
%|    #    |   FILE       : mshNdvolume.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Element volume of dimension n                 |
%|  `---'  |                                                              |
%+========================================================================+
 
% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshNdvolume.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    
    % Volume
    V = sqrt(sum(E1.^2,2));
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
    E3 = cross(E1,E2,2);
    
    % Volume
    V = 0.5*sqrt(sum(E3.^2,2));
    
% Tetrahedral mesh
elseif (size(mesh.elt,2) == 4)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
    E3 = mesh.vtx(mesh.elt(:,4),:) - mesh.vtx(mesh.elt(:,1),:);

    % Volume
    V = 1/6 .* abs( ...
          E1(:,1) .* (E2(:,2).*E3(:,3) - E2(:,3).*E3(:,2)) ...
        - E2(:,1) .* (E1(:,2).*E3(:,3) - E1(:,3).*E3(:,2)) ...
        + E3(:,1) .* (E1(:,2).*E2(:,3) - E1(:,3).*E2(:,2)) );
    
% Unknown type
else
    error('mshNdvolume.m : unavailable case')
end
end
