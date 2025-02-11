function [mesh,elt2edg] = mshEdge(mesh)
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
%|    #    |   FILE       : mshEdge.m                                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Submesh to edge mesh                          |
%|  `---'  |                                                              |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshEdge.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    edg2vtx = mesh.elt;
    elt2edg = (1:size(mesh.elt,1))';
    col     = mesh.col;
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    % All edges
    edg2vtx = [ mesh.elt(:,[2 3]) ; mesh.elt(:,[3 1]) ; mesh.elt(:,[1 2]) ];
    
    % Edges unicity
    tmp           = sort(edg2vtx,2);
    [~,I,elt2edg] = unique(tmp,'rows');
    
    % Final sorted edges 
    edg2vtx = tmp(I,:);
    elt2edg = reshape(elt2edg,size(mesh.elt,1),3);
    
    % Colours
    col             = zeros(size(edg2vtx,1),1);
    tmp             = mesh.col * ones(1,3);
    col(elt2edg(:)) = tmp(:);
    
    % Identify boundary
    bnd      = ( accumarray(elt2edg(:),1).*col ~= accumarray(elt2edg(:),tmp(:)) );
    col(bnd) = -1;    
    
% Tetrahedral mesh
elseif (size(mesh.elt,2) == 4)
    % All edges
    edg2vtx = [ mesh.elt(:,[1 2]) ; mesh.elt(:,[1 3]) ; mesh.elt(:,[1 4]) ; 
        mesh.elt(:,[2 3]) ; mesh.elt(:,[2 4]) ; mesh.elt(:,[3 4]) ]; 
        
    % Edges uniqueness
    tmp           = sort(edg2vtx,2);
    [~,I,elt2edg] = unique(tmp,'rows');
    
    % Final Edges
    edg2vtx = tmp(I,:);
    elt2edg = reshape(elt2edg,size(mesh.elt,1),6);
    
    % Colours
    col             = zeros(size(edg2vtx,1),1);
    tmp             = mesh.col * ones(1,6);
    col(elt2edg(:)) = tmp(:);
    
    % Identify boundary
    bnd      = ( accumarray(elt2edg(:),1).*col ~= accumarray(elt2edg(:),tmp(:)) );
    col(bnd) = -1;
    
% Unknown type
else
    error('mshEdge.m : unavailable case')
end

% Mesh format
mesh = msh(mesh.vtx,edg2vtx,col);
end
