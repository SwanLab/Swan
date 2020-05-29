function [X,elt2dof] = femDof(fe)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : femDof.m                                      |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Define degrees of freedom                     |
%|  `---'  |                                                              |
%+========================================================================+

% Lagrange order 0, constant by element
if strcmp(fe.typ,'P0')
    X       = fe.msh.ctr;
    elt2dof = (1:size(fe.msh.elt,1))';
    
% Lagrange order 1, piecewise linear by element
elseif strcmp(fe.typ,'P1')
    X       = fe.msh.vtx;
    elt2dof = fe.msh.elt;
    
% Lagrange order 2, piecewise quadratic by element
elseif strcmp(fe.typ,'P2')
    [edge,elt2edg] = fe.msh.edg;
    X              = [fe.msh.vtx ; edge.ctr];
    elt2dof        = [fe.msh.elt , elt2edg + size(fe.msh.vtx,1)];
    
% NED linear elements
elseif strcmp(fe.typ,'NED')
    [edge,elt2dof] = fe.msh.edg;
    X              = edge.ctr;
    
% RWG linear elements for surfacic mesh
elseif strcmp(fe.typ,'RWG')
    % Particular mesh (dof are vertices)
    if (size(fe.msh.elt,2) == 1)
        X       = fe.msh.vtx;
        elt2dof = fe.msh.elt;
        
    % Triangular mesh
    elseif (size(fe.msh.elt,2) == 3)
        [edge,elt2dof] = fe.msh.edg;
        X              = edge.ctr;
        
    % Tetrahedral mesh
    elseif size(fe.msh.elt,2)==4
        [face,elt2dof] = fe.msh.fce;
        X              = face.ctr;
    end
    
% Others
else
    error('femDof.m : unavailable case')
end
end
