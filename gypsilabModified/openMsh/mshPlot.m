function mshPlot(mesh,V)
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
%|    #    |   FILE       : mshPlot.m                                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Plot mesh and data                            |
%|  `---'  |                                                              |
%+========================================================================+

% Patch 
H = patch('Faces',mesh.elt, 'Vertices',mesh.vtx,'Marker','o',...
    'MarkerFaceColor','none','EdgeColor','none','FaceColor','none');

% For single numerical values
if (numel(V) == 1) && ~ischar(V)
    V = V * ones(size(mesh.vtx,1),1);
end

% Particles mesh
if (size(mesh.elt,2) == 1)
    if isempty(V)
        mlt = accumarray(mesh.elt,1);
        col = accumarray(mesh.elt,mesh.col)./mlt;
        set(H,'MarkerFaceColor','flat','FaceVertexCData',col);
    elseif ischar(V)
        set(H,'MarkerFaceColor',V);
    else
        set(H,'MarkerFaceColor','flat','FaceVertexCData',V);
    end   

% Edge mesh
elseif (size(mesh.elt,2) == 2)
    if isempty(V)
        Nelt = size(mesh.elt,1);
        vtx  = [mesh.vtx(mesh.elt(:,1),:) ; mesh.vtx(mesh.elt(:,2),:)];
        elt  = [(1:Nelt)' (Nelt+1:2*Nelt)'];
        col  = [mesh.col ; mesh.col];
        delete(H);
        patch('Faces',elt, 'Vertices',vtx,'Marker','o',...
            'MarkerFaceColor','k','EdgeColor','flat','FaceColor','none',...
            'FaceVertexCData',col);
    elseif ischar(V)
        set(H,'MarkerFaceColor','k','EdgeColor',V);
    else
        set(H,'Marker','none','EdgeColor','interp','FaceVertexCData',V);
    end
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    if isempty(V)
        set(H,'Marker','.','MarkerFaceColor','k','EdgeColor','k','FaceColor','flat','CData',mesh.col)
    elseif ischar(V)
        set(H,'Marker','.','MarkerFaceColor','k','EdgeColor','k','FaceColor',V)
    else
        set(H,'Marker','none','EdgeColor','none','FaceColor','interp','FaceVertexCData',V)
    end    
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    mshPlot(mesh.fce,V);
    
% Unknown type    
else
    error('mshPlot.m : unavailable case')
end
end
