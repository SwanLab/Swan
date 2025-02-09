function [X,P] = femUnk(fe)
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
%|    #    |   FILE       : femUnk.m                                      |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Unknowns and reduction matrix for constrained |
%|  `---'  |                finite elements                               |
%+========================================================================+

% Initialization with dofs
X  = fe.dof;
Nx = size(X,1);
P  = speye(Nx);

% Dirichlet
if ~isempty(fe.dir)
    % Finite element for dirichlet mesh
    fed = fem(fe.dir,fe.typ);
    
    % Dof indices
    I = find(~ismember(X,fed.dof,'rows'));
    
    % Reduction matrix
    M = speye(Nx);
    M = M(:,I);
    
    % Update
    P  = P * M;
    X  = X(I,:);
    Nx = length(I);
end

% Junction
if ~isempty(fe.jct)
    % Number of junction
    Njct = length(fe.jct)/2;
    
    % Initialization
    feJct = cell(1,Njct);
    Ijct  = cell(1,Njct);
    
    % Finite element for junction meshes
    for i = 1:Njct
        if isa(fe.jct{2*i-1},'msh') && isnumeric(fe.jct{2*i})
            feJct{i} = fem(fe.jct{2*i-1},fe.typ);
        else
            error('femUnknown.m : unavailable case');
        end
    end
    
    % Valid dof indices for junction meshes    
    for i = 1:Njct
        Xi         = feJct{i}.dof;
        Ijct{i}    = zeros(size(Xi,1),1);
        [~,I,J]    = intersect(Xi,X,'rows','stable');
        Ijct{i}(I) = J; 
    end
    Ijct = cell2mat(Ijct);
     
    % Extract dirichlet condition and multiple indices
    bool = (Ijct(:,Njct) == 0);
    for i = 1:Njct-1
        bool = bool + (Ijct(:,Njct) == Ijct(:,i));
    end
    Ijct = Ijct(~bool,:);
    
    % Final dof after extraction of the last indices
    I = setdiff((1:Nx)',Ijct(:,end));

    % Linear relation matrix
    M = speye(Nx);
    for  i = 1:Njct
        for j = setdiff(1:Njct,i)            
            M = M + sparse(Ijct(:,i),Ijct(:,j),-fe.jct{2*j}/fe.jct{2*i},Nx,Nx);
        end
    end    
    
    % Reduction
    M = M(:,I);
    
    % Update
    P  = P * M;
    X  = X(I,:);
%     Nx = length(I);
end
end
