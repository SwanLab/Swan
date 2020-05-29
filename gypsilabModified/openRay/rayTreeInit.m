function tr = rayTreeInit(ray)
%+========================================================================+
%|                                                                        |
%|           OPENRAY - LIBRARY FOR TRI-DIMENSIONAL RAY TRACING            |
%|           openRay is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : rayTreeInit.m                                 |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Prepare tree for ray propagation              |
%|  `---'  |                                                              |
%+========================================================================+

% Tree
tr = tree(ray.msh,'binary');

% Loop at each step of the tree
for n = 1:length(tr)
    % Leaf size
    Nleaf = length(tr{n}.ind);
    
    % Initialize children ray indices
    tr{n}.rad = zeros(Nleaf,1);
    
    % Update data
    for i = 1:Nleaf
        % Element indices
        ielt = tr{n}.ind{i};

        % Local vertices
        elt = ray.msh.elt(ielt,:);
        vtx = ray.msh.vtx(elt(:),:);
        
        % Update center
        M              = 0.5 * (min(vtx,[],1) + max(vtx,[],1));
        tr{n}.ctr(i,:) = M;
        
        % Update spherical radius
        tr{n}.rad(i) = sqrt( max( (vtx(:,1)-M(1)).^2 + ...
            (vtx(:,2)-M(2)).^2 + ...
            (vtx(:,3)-M(3)).^2 ) );
    end
end
end
