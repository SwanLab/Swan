function mesh = mshDisk(N,rad)
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
%|    #    |   FILE       : mshDisk.m                                     |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Build uniform mesh for a disk                 |
%|  `---'  |                                                              |
%+========================================================================+

% Radial discretisation
dr = sqrt(pi*rad^2/N);
dr = rad/ceil(rad/dr);
r  = dr:dr:rad;

% Angular uniform discretization
rho = cell(length(r),1); theta = rho;
for ir = 1:length(r)
    dtheta = dr/r(ir);
    dtheta = 2*pi/ceil(2*pi/dtheta);    
    theta{ir} = (0:dtheta:2*pi-dtheta)';
    rho{ir}   = r(ir)*ones(length(theta{ir}),1);
end

% Carthesian coordinates
[x,y] = pol2cart(cell2mat(theta),cell2mat(rho));
X     = [0 0 ; x y];

% Unicity test
tmp = unique(X,'rows','stable');
if (max(abs(X-tmp)) > 1e-12)
    error('mshDisk : non unicity of vertices')
end
   
% Delaunay triangulation
DT = delaunayTriangulation(X(:,1),X(:,2));

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(vtx,elt);
end
