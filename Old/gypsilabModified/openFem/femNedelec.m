function M = femNedelec(fe,domain)
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
%|    #    |   FILE       : femNedelec.m                                  |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Nedelec finite element matrix                 |
%|  `---'  |                                                              |
%+========================================================================+

% Gaussian quadrature
Xgss             = domReference(domain);
[Xqud,~,elt2qud] = domain.qud;

% Edges
[mshedg,elt2edg] = fe.msh.edg;

% Intersect domain and finite element meshes
if isequal(fe.msh,domain.msh)
    mesh = fe.msh;
    Ife  = (1:size(fe.msh.elt,1))';
    Idom = Ife;
    XintQud = Xqud;
else
    [mesh,Ife,Idom] = intersect(fe.msh,domain.msh);
    domainInt       = dom(mesh,domain.gss);
    XintQud         = domainInt.qud;
    
end

% Dimensions
Nelt = length(Ife);
Nqud = size(Xqud,1);
Ngss = size(Xgss,1);
% Dof are the edges of the mesh
Ndof = size(mshedg.elt,1);
Nbas = size(elt2edg,2);
% Volumes
Ndv  = mesh.ndv;

% Numbering dof to quadrature
idx = zeros(Nelt,Nbas,Ngss);
jdx = zeros(Nelt,Nbas,Ngss);
for i = 1:Nbas
    for j = 1:Ngss
        idx(:,i,j) = elt2qud(Idom,j);
        jdx(:,i,j) = elt2edg(Ife,i);
    end
end
        
%%% FINITE ELEMENT MATRIX AND GRADIENT
if size(fe.msh.elt,2)==3 % Triangular elements
    
    % Vectorial operator nx[PSI]
    if strcmp(fe.opr,'nx[psi]')
        M = cell(1,3);
        for n = 1:3
            % Initialization
            val = zeros(Nelt,Nbas,Ngss);
            
            % For each basis function
            for i = 1:Nbas
                % Indices of neighbors
                ip1 = mod(i,3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the edge
                flux = (2*(mesh.elt(:,ip1) < mesh.elt(:,ip2))-1)./(2*Ndv);
                
                % For each integration point
                for j = 1:Ngss
                    val(:,i,j) = flux.*(XintQud(j:Ngss:Nqud,n) - mesh.vtx(mesh.elt(:,i),n));
                end
            end
            
            % Operator
            M{n} = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
        end
        
        % Scalar operator div[PSI]
    elseif strcmp(fe.opr,'curl[psi]')
        % Initialization
        val = zeros(Nelt,Nbas,Ngss);
        
        % For each basis function
        for i = 1:Nbas
            % Indices of neighbors
            ip1 = mod(i,3)+1;
            ip2 = mod(ip1,3)+1;
            
            % Flux through the edge
            flux = (2*(mesh.elt(:,ip1) < mesh.elt(:,ip2))-1)./Ndv;
            
            % For each integration point
            for j = 1:Ngss
                val(:,i,j) = flux;
            end
        end
        
        % Operator
        M = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% NORMAL x DQM
    elseif strcmp(fe.opr,'[psi]')
        % Finite element
        tmp     = fem(fe.msh,fe.typ);
        tmp.opr = 'nx[psi]';
        uqm     = tmp.uqm(domain);
        
        % Normals
        nrm  = domain.qudNrm;
        m    = size(nrm,1);
        N{1} = spdiags(nrm(:,1),0,m,m);
        N{2} = spdiags(nrm(:,2),0,m,m);
        N{3} = spdiags(nrm(:,3),0,m,m);
        
        % Cross product
        M = cell(1,3);
        for i = 1:3
            ip1  = mod(i,3) + 1;
            ip2  = mod(ip1,3) + 1;
            M{i} = - (N{ip1} * uqm{ip2} - N{ip2} * uqm{ip1});
        end
        
    else
        error('femNedelec.m : unavailable case')
    end
elseif size(fe.msh.elt,2)==4 % Tetrahedral elements
    
    M = cell(1,3);
    for n = 1:3
        % Initialisation
        val = zeros(Nelt,Nbas,Ngss);
        
        % Indices des coordonnees
        np1 = mod(n,3)+1;
        np2 = mod(np1,3)+1;
        
        % For each base function
        for i = 1:Nbas
            % Opposite edge
            ibar = 7-i;
            % 1st vertex of the edge
            xi = mesh.vtx(mshedg.elt(elt2edg(:,i),1),:);
            % 1st vertex of opposite edge
            xibar = mesh.vtx(mshedg.elt(elt2edg(:,ibar),1),:);
            % Current edge
            ai = mesh.vtx(mshedg.elt(elt2edg(:,i),2),:) - xi;
            % opposite edge
            aibar = mesh.vtx(mshedg.elt(elt2edg(:,ibar),2),:) - xibar;
            % Correcting coefficient
            bi = xi - xibar;
            alphai = 1./(ai(:,1).*aibar(:,2).*bi(:,3) - ai(:,1).*aibar(:,3).*bi(:,2) ...
                + ai(:,2).*aibar(:,3).*bi(:,1) - ai(:,2).*aibar(:,1).*bi(:,3) ...
                + ai(:,3).*aibar(:,1).*bi(:,2) - ai(:,3).*aibar(:,2).*bi(:,1));
            if strcmp(fe.opr,'[psi]')
                % Loop over integration points
                for j = 1:Ngss
                    x = XintQud(j:Ngss:Nqud,:);
                    val(:,i,j) = ...
                        alphai.*(aibar(:,np1).*(x(:,np2) - xibar(:,np2)) - aibar(:,np2).*(x(:,np1) - xibar(:,np1)));
                end
            elseif strcmp(fe.opr,'curl[psi]')
                % Loop over integration points
                for j = 1:Ngss
                    val(:,i,j) = 2 * alphai .* aibar(:,n);
                end
            else
                error('femNedelec.m : unavailable case')
            end
        end
        M{n} = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
    end
else
    error('femNedelec.m : unavailable case')
end
end
