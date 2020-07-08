function [X,W,elt2qud] = domQuadrature(domain)
%%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : domQuadrature.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Discrete quadrature for edges, triangle       |
%|  `---'  |                and tetrahedre meshes                         |
%+========================================================================+
    
% Reference quadrature
[x,w] = domReference(domain);

% Initialization
X       = zeros(size(x,1)*size(domain.msh.elt,1),size(domain.msh.vtx,2));
elt2qud = zeros(size(domain.msh.elt,1),size(x,1));

% For each quadrature point
for j = 1:size(x,1)
    % Indice
    idx          = (j:size(x,1):size(x,1)*size(domain.msh.elt,1))';
    elt2qud(:,j) = idx;
    
    % Particle mesh
    if (size(domain.msh.elt,2) == 1)
        error('domQuadrature.m : unavailable case')
            
    % Edge mesh
    elseif (size(domain.msh.elt,2) == 2)
        X(idx,:) = (1-x(j,1)) * domain.msh.vtx(domain.msh.elt(:,1),:) ...
            + x(j,1) * domain.msh.vtx(domain.msh.elt(:,2),:) ;

    % Triangular mesh
    elseif (size(domain.msh.elt,2) == 3)
        X(idx,:) = (1-x(j,1)-x(j,2)) * domain.msh.vtx(domain.msh.elt(:,1),:) ...
            + x(j,1) * domain.msh.vtx(domain.msh.elt(:,2),:) ...
            + x(j,2) * domain.msh.vtx(domain.msh.elt(:,3),:);
        
    % Tetrahedral mesh
    elseif (size(domain.msh.elt,2) == 4)
        X(idx,:) = (1-x(j,1)-x(j,2)-x(j,3)) * domain.msh.vtx(domain.msh.elt(:,1),:) ...
            + x(j,1) * domain.msh.vtx(domain.msh.elt(:,2),:) ...
            + x(j,2) * domain.msh.vtx(domain.msh.elt(:,3),:) ...
            + x(j,3) * domain.msh.vtx(domain.msh.elt(:,4),:);
        
    % Unknown type
    else
        error('domQuadrature.m : unavailable case')
    end
end

% Quadrature weight
W = kron(domain.msh.ndv,w);
end
