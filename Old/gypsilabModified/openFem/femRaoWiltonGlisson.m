function M = femRaoWiltonGlisson(fe,domain)
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
%|    #    |   FILE       : femRaoWiltonGlisson.m                         |
%|    #    |   VERSION    : 0.60                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : RWG finite element matrix                     |
%|  `---'  |                                                              |
%+========================================================================+

%%% FINITE ELEMENT MATRIX AND DIVERGENCE
if strcmp(fe.opr,'[psi]') || strcmp(fe.opr,'div[psi]')
    % Gaussian quadrature
    Xgss             = domReference(domain);
    [Xqud,~,elt2qud] = domain.qud;
    
    % Degrees of freedom
    [Xdof,elt2dof] = fe.dof;
    
    % Intersect domain and finite element meshes
    if isequal(fe.msh,domain.msh)
        mesh    = fe.msh;
        Ife     = (1:size(fe.msh.elt,1))';
        Idom    = Ife;
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
    Ndof = size(Xdof,1);
    Nbas = size(elt2dof,2);
    Ndv  = mesh.ndv;
    
    % Numbering dof to quadrature
    idx = zeros(Nelt,Nbas,Ngss);
    jdx = zeros(Nelt,Nbas,Ngss);
    for i = 1:Nbas
        for j = 1:Ngss
            idx(:,i,j) = elt2qud(Idom,j);
            jdx(:,i,j) = elt2dof(Ife,i);
        end
    end
end

%%% FINITE ELEMENT MATRIX AND GRADIENT
if size(fe.msh.elt,2)==3 % Triangular elements
    % Vectorial operator [PSI]
    if strcmp(fe.opr,'[psi]')
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
                    val(:,i,j) = flux.*(XintQud(j:Ngss:end,n) - mesh.vtx(mesh.elt(:,i),n));
                end
            end
            
            % Operator
            M{n} = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
        end
        
        % Scalar operator div[PSI]
    elseif strcmp(fe.opr,'div[psi]')
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
        %%% RESTRICTION TO Ith COMPONENT
    elseif strcmp(fe.opr(1:5),'[psi]')
        % Finite element
        tmp     = fem(fe.msh,fe.typ);
        tmp.opr = '[psi]';
        uqm     = tmp.uqm(domain);
        
        % Restriction
        ind     = str2double(fe.opr(end));
        M       = uqm{ind};

        
        %%%% NORMAL x DQM
    elseif strcmp(fe.opr,'nx[psi]')
        % Finite element
        tmp     = fem(fe.msh,fe.typ);
        tmp.opr = '[psi]';
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
            M{i} = N{ip1} * uqm{ip2} - N{ip2} * uqm{ip1};
        end
        
        
        %%%% NORMAL x DQM (j)
    elseif strcmp(fe.opr(1:end-1),'nx[psi]')
        % Component
        j = str2double(fe.opr(end));
        
        % Finite element
        tmp     = fem(fe.msh,fe.typ);
        tmp.opr = '[psi]';
        uqm     = tmp.uqm(domain);
        
        % Normals
        N = domain.qudNrm;
        m = size(N,1);
        
        % Cross component
        jp1 = mod(j,3) + 1;
        jp2 = mod(jp1,3) + 1;
        
        % Cross product
        M = spdiags(N(:,jp1),0,m,m) * uqm{jp2} - ...
            spdiags(N(:,jp2),0,m,m) * uqm{jp1};
        
    else
        error('femRaoWiltonGlisson.m : unavailable case')
    end
elseif size(fe.msh.elt,2)==4 % Tetrahedral elements
    % Vectorial operator [PSI]
    if strcmp(fe.opr,'[psi]')
        M = cell(1,3);
        for n = 1:3
            % Initialization
            val = zeros(Nelt,Nbas,Ngss);
            
            % For each basis function
            for i = 1:Nbas
                % Indices of neighbors
                ip1 = mod(i,4)+1;
                ip2 = mod(ip1,4)+1;
                ip3 = mod(ip2,4)+1;
                
                % Flux through the face
                flux = (2*(mesh.elt(:,ip1) < mesh.elt(:,ip2))-1) ...
                    .*(2*(mesh.elt(:,ip2) < mesh.elt(:,ip3))-1) ...
                    .*(2*(mesh.elt(:,ip1) < mesh.elt(:,ip3))-1) ...
                    ./(3*Ndv);
                if mod(i,2)==0
                    flux = - flux;
                end
                
                % For each integration point
                for j = 1:Ngss
                    val(:,i,j) = flux.*(XintQud(j:Ngss:end,n) - mesh.vtx(mesh.elt(:,i),n));
                end
            end
            
            % Operator
            M{n} = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
        end
        
    % Scalar operator div[PSI]
    elseif strcmp(fe.opr,'div[psi]')
        % Initialization
        val = zeros(Nelt,Nbas,Ngss);
        
        % For each basis function
        for i = 1:Nbas
            % Indices of neighbors
            ip1 = mod(i,4)+1;
            ip2 = mod(ip1,4)+1;
            ip3 = mod(ip2,4)+1;
            
            % Flux through the edge
            flux = (2*(mesh.elt(:,ip1) < mesh.elt(:,ip2))-1) ...
                .*(2*(mesh.elt(:,ip2) < mesh.elt(:,ip3))-1) ...
                .*(2*(mesh.elt(:,ip1) < mesh.elt(:,ip3))-1) ...
                ./Ndv;
            if mod(i,2)==0
                flux = - flux;
            end
            
            % For each integration point
            for j = 1:Ngss
                val(:,i,j) = flux;
            end
        end
        
        % Operator
        M = sparse(idx(:),jdx(:),val(:),Nqud,Ndof);
    else
        error('femRaoWiltonGlisson.m : unavailable case')
    end
end
end
