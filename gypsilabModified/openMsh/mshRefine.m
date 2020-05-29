function mesh = mshRefine(mesh,ord)
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
%|    #    |   FILE       : mshRefine.m                                   |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Refine triangular mesh with various order of  |
%|  `---'  |                midpoint algorithm                            |
%+========================================================================+

% Fixed order for all element
if (isnumeric(ord) || islogical(ord)) && (length(ord)==length(mesh))
    % Recursive loop
    while (sum(ord) > 0)
        % Midpoint refinement
        [mesh,Ir] = midpoint(mesh,find(ord));
        
        % Update order
        ord = ord(Ir) - 1;
        ord = ord .* (ord > 0);
    end

% Fixed edge length    
elseif isnumeric(ord) && (length(ord)==1)
    % Select element from edge mesh
    [edg,Itri] = mesh.edg;
    Ilgt       = (ndv(edg) > ord);
    I          = sum(Ilgt(Itri),2);
    
    % Recursive loop
    while (sum(I) > 0)
        % Midpoint refinement
        mesh = midpoint(mesh,find(I));
        
        % Select element from edge mesh
        [edg,Itri] = mesh.edg;
        Ilgt       = ndv(edg) > ord;
        I          = sum(Ilgt(Itri),2);
    end

% Function order
elseif isa(ord,'function_handle')
    % Initialize indices
    I = ord(mesh.ctr);
    n = 1;

    % Recursive loop
    while (sum(I) > 0)
        % Midpoint refinement
        mesh = midpoint(mesh,find(I));
        
        % Update order
        I = ord(mesh.ctr) - n;
        I = I .* (I > 0);
        n = n + 1;
    end
    
% Unknown type
else
    error('mshRefine.m : unknown case')
end
end
