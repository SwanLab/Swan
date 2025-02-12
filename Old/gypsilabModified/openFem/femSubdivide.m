function [Isub,domSub,femSub] = femSubdivide(dom,fem,Nsub,fig)
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
%|    #    |   FILE       : femSubdivide.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Subdivide finite element space for domain     |
%|  `---'  |                decomposition method                          |
%|         |                ! ONLY WORKS FOR NON DIRICHLET FEM !          |
%+========================================================================+

% Number of subdomains is power of 2
if (Nsub ~= 2^floor(log2(Nsub)))
    error('femSubdivide.m : number of subdomain is not a power of 2.')
end

% Mesh and dof associated to finite element space
mesh          = fem.msh;
[dof,elt2dof] = fem.dof;                                   % to be improved

% Dof subdivision with binary tree
Nleaf  = ceil(length(fem)/Nsub);
bitree = tree(msh(dof),'binary',Nleaf);
Isub   = bitree{end}.ind;

% Security
if (length(Isub) ~= Nsub)
    error('femSubdivide.m : unavailable case')
end
if (norm(sort(cell2mat(Isub)) - (1:length(fem))','inf') > 1e-12)
    error('femSubdivide.m : unavailable case')    
end

% Output initialization
domSub = cell(Nsub,1);
femSub = cell(Nsub,1);

% For each subindices
for i = 1:Nsub
    % Indices of valid dof and extended mesh
    I       = ismember(elt2dof,Isub{i});
    meshSub = mesh.sub(sum(I,2)>0);
    
    % Quadrature
    domSub{i}     = dom;
    domSub{i}.msh = meshSub;
    
    % Finite element
    femSub{i}     = fem;
    femSub{i}.msh = meshSub;
    
    % Dirichlet condition for boundary dof
    dir = setdiff( msh(femSub{i}.dof) , msh(dof(Isub{i},:)) );
    if (size(dir,1)>0)
        femSub{i} = dirichlet(femSub{i},dir);
    end
    
    % Security
    if (length(femSub{i}) ~= length(Isub{i}))
        error('femSubdivide.m : unavailable case');
    end

    % Graphical representation
    if fig
        tmp = mesh.sub(sum(I,2)==size(elt2dof,2));
        figure(fig)
        hold on
        plot(tmp,'b')
        plot(tmp.bnd,'r')
        plot(setdiff(meshSub,tmp),'w')        
        axis equal
    end
end
end
