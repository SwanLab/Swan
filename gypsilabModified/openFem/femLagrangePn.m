function M = femLagrangePn(fe,domain)
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
%|    #    |   FILE       : femLagrangePn.m                               |
%|    #    |   VERSION    : 0.55                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.05.2019                                    |
%| ( === ) |   SYNOPSIS   : Lagrange finite element matrix (P0,P1,P2)     |
%|  `---'  |                                                              |
%+========================================================================+

%%% FINITE ELEMENT MATRIX AND GRADIENT
if strcmp(fe.opr,'[psi]') || strcmp(fe.opr,'grad[psi]')
    % Gaussian quadrature
    Xgss             = domReference(domain);
    [Xqud,~,elt2qud] = domain.qud;
    
    % Degrees of freedom
    [Xdof,elt2dof] = fe.dof;
    
    % Intersect domain and finite element meshes
    if isequal(fe.msh,domain.msh)
        mesh = fe.msh;
        Ife  = (1:size(fe.msh.elt,1))';
        Idom = Ife;
    else
        [mesh,Ife,Idom] = intersect(fe.msh,domain.msh);
    end
    
    % Dimensions
    Nelt = length(Ife);
    Nqud = size(Xqud,1);
    Ngss = size(Xgss,1);
    Ndof = size(Xdof,1);
    Nbas = size(elt2dof,2);
    dim  = size(fe.msh.elt,2);    
        
    % Numbering dof to quadrature
    idx = zeros(Nelt,Nbas,Ngss);
    jdx = zeros(Nelt,Nbas,Ngss);
    for i = 1:Nbas
        for j = 1:Ngss
            idx(:,i,j) = elt2qud(Idom,j);
            jdx(:,i,j) = elt2dof(Ife,i);
        end
    end
    
    % Basis fct constant per elements
    if strcmp(fe.typ,'P0')
        % Edge mesh
        if (dim == 2)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);

        % Triangular mesh
        elseif (dim == 3)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);
            dyF = zeros(1,Ngss);
            
        % Tetrahedron mesh
        elseif (dim == 4)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);
            dyF = zeros(1,Ngss);
            dzF = zeros(1,Ngss);
        
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    % Basis fct piecewise linear per element
    elseif strcmp(fe.typ,'P1')
        % Initialization
        F   = zeros(dim,Ngss);
        dxF = zeros(dim,Ngss);
        dyF = zeros(dim,Ngss);
        dzF = zeros(dim,Ngss) ;       
        
        % Edge mesh
        if (dim == 2)
            F(1,:) = 1 - Xgss(:,1);
            F(2,:) = Xgss(:,1);
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;

        % Triangular mesh
        elseif (dim == 3)
            F(1,:) = 1 - Xgss(:,1) - Xgss(:,2);
            F(2,:) = Xgss(:,1);
            F(3,:) = Xgss(:,2);            
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;
            
            dyF(1,:) = - 1;
            dyF(3,:) = 1;
            
        % Tetrahedron mesh
        elseif (dim == 4)
            F(1,:) = 1 - Xgss(:,1) - Xgss(:,2) - Xgss(:,3);
            F(2,:) = Xgss(:,1);
            F(3,:) = Xgss(:,2);
            F(4,:) = Xgss(:,3);
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;
                        
            dyF(1,:) = - 1;
            dyF(3,:) = 1;
            
            dzF(1,:) = - 1;
            dzF(4,:) = 1;
            
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    % Lagrange order 2, piecewise quadratic by element     
    elseif strcmp(fe.typ,'P2')
        % Edge mesh
        if (size(mesh.elt,2) == 2)
            % Initialization
            F   = zeros(3,Ngss);
            dxF = zeros(3,Ngss);
            dyF = zeros(3,Ngss);
            dzF = zeros(3,Ngss);
            X = Xgss(:,1);
            
            F(1,:) = (1 - X).*(1 - 2*X);
            F(2,:) = X.*(2*X - 1);
            F(3,:) = X.*(1 - X);
            
            dxF(1,:) = 4*X - 3;
            dxF(2,:) = 4*X - 1;
            dxF(3,:) = 2*X - 1;
            
        % Triangular mesh
        elseif (size(mesh.elt,2) == 3)
            % Initialization
            F   = zeros(6,Ngss);
            dxF = zeros(6,Ngss);
            dyF = zeros(6,Ngss);
            dzF = zeros(6,Ngss);
            X = Xgss(:,1);
            Y = Xgss(:,2);
            
            F(1,:) = (1 - X - Y).*(1 - 2*X - 2*Y);
            F(2,:) = X.*(2*X - 1);
            F(3,:) = Y.*(2*Y - 1);
            F(4,:) = 4*X.*Y;
            F(5,:) = 4*Y.*(1 - X - Y);
            F(6,:) = 4*X.*(1 - X - Y);
            
            dxF(1,:) = -3 + 4*(X + Y);
            dxF(2,:) = 4*X - 1;
            dxF(4,:) = 4*Y;
            dxF(5,:) = -4*Y;
            dxF(6,:) = 4*(1 - 2*X - Y);
            
            dyF(1,:) = -3 + 4*(X + Y);
            dyF(3,:) = 4*Y - 1;
            dyF(4,:) = 4*X;
            dyF(5,:) = 4*(1 - X - 2*Y);
            dyF(6,:) = -4*X;
            
        % Tetrahedron mesh
        elseif (size(mesh.elt,2) == 4)
            % Initialization
            F   = zeros(10,Ngss);
            dxF = zeros(10,Ngss);
            dyF = zeros(10,Ngss);
            dzF = zeros(10,Ngss);
            X = Xgss(:,1);
            Y = Xgss(:,2);
            Z = Xgss(:,3);
            
            F(1,:) = (1 - X - Y - Z).*(1 - 2*X - 2*Y - 2*Z);
            F(2,:) = X.*(2*X - 1);
            F(3,:) = Y.*(2*Y - 1);
            F(4,:) = Z.*(2*Z - 1);
            F(5,:) = 4*X.*(1 - X - Y - Z);
            F(6,:) = 4*X.*Y;
            F(7,:) = 4*Y.*(1 - X - Y - Z);
            F(9,:) = 4*Z.*(1 - X - Y - Z);
            F(8,:) = 4*X.*Z;
            F(10,:) = 4*Y.*Z;
             
            dxF(1,:) = -3 + 4*(X + Y + Z);
            dxF(2,:) = 4*X - 1;
            dxF(5,:) = 4*(1 - 2*X - Y - Z);
            dxF(6,:) = 4*Y;
            dxF(7,:) = -4*Y;
            dxF(9,:) = -4*Z;
            dxF(8,:) = 4*Z;
                        
            dyF(1,:) = -3 + 4*(X + Y + Z);
            dyF(3,:) = 4*Y - 1;
            dyF(5,:) = -4*X;
            dyF(6,:) = 4*X;
            dyF(7,:) = 4*(1 - X - 2*Y - Z);
            dyF(9,:) = -4*Z;
            dyF(10,:) = 4*Z;

            dzF(1,:) = -3 + 4*(X + Y + Z);
            dzF(4,:) = 4*Z - 1;
            dzF(5,:) = -4*X;
            dzF(7,:) = -4*Y;
            dzF(9,:) = 4*(1 - X - Y - 2*Z);
            dzF(8,:) = 4*X;
            dzF(10,:) = 4*Y;
            
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    else
        error('femLagrangePn.m: unavailable case');
    end
    
    % Finite element matrix
    if strcmp(fe.opr,'[psi]')
        bas = zeros(Nelt,Nbas,Ngss);
        for i = 1:Nbas
            for j = 1:Ngss
                bas(:,i,j) = F(i,j);
            end
        end
        M = sparse(idx(:),jdx(:),bas(:),Nqud,Ndof);
        
    % Gradient of Finite element matrix
    elseif strcmp(fe.opr,'grad[psi]')
        % Edge mesh
        if (dim == 2)
            % Vector basis
            E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            
            % Inverse Gramm matrix
            Dx1 = 1./(E1(:,1).^2 + E1(:,2).^2 + E1(:,3).^2);
            
            % Gradient projection to integration points
            dbas = cell(1,3);
            for n = 1:3
                dbas{n} = zeros(Nelt,Nbas,Ngss);
                DCVx    = Dx1.*E1(:,n);
                for i = 1:Nbas
                    for j = 1:Ngss
                        dbas{n}(:,i,j) = DCVx(:)*dxF(i,j);
                    end
                end
            end
            
        % Triangular mesh
        elseif (dim == 3)
            % Vector basis
            E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
            
            % Gramm matrix coefficients
            a = E1(:,1).^2 + E1(:,2).^2 + E1(:,3).^2;
            b = E1(:,1).*E2(:,1) + E1(:,2).*E2(:,2) + E1(:,3).*E2(:,3);
            c = b;
            d = E2(:,1).^2 + E2(:,2).^2 + E2(:,3).^2;
            
            % Determinant
            detG = a.*d - b.*c;
            
            % Inverse by co-factor
            Dx1 = d  ./ detG;
            Dx2 = -b ./ detG;
            Dy1 = -c ./ detG;
            Dy2 = a  ./ detG;
            
            % Gradient projection to integration points
            dbas = cell(1,3);
            for n = 1:3
                dbas{n} = zeros(Nelt,Nbas,Ngss);
                DCVx    = Dx1.*E1(:,n) + Dx2.*E2(:,n);
                DCVy    = Dy1.*E1(:,n) + Dy2.*E2(:,n);
                for i = 1:Nbas
                    for j = 1:Ngss
                        dbas{n}(:,i,j) = DCVx(:)*dxF(i,j) + DCVy(:)*dyF(i,j); 
                    end
                end
            end
            
        % Tetrahedron mesh
        elseif (dim == 4)            
            % Vector basis
            E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
            E3 = mesh.vtx(mesh.elt(:,4),:) - mesh.vtx(mesh.elt(:,1),:);
            
            % Gramm matrix
            a = E1(:,1).^2 + E1(:,2).^2 + E1(:,3).^2;
            b = E1(:,1).*E2(:,1) + E1(:,2).*E2(:,2) + E1(:,3).*E2(:,3);
            c = E1(:,1).*E3(:,1) + E1(:,2).*E3(:,2) + E1(:,3).*E3(:,3);
            d = b;
            e = E2(:,1).^2 + E2(:,2).^2 + E2(:,3).^2;
            f = E2(:,1).*E3(:,1) + E2(:,2).*E3(:,2) + E2(:,3).*E3(:,3);
            g = c;
            h = f;
            i = E3(:,1).^2 + E3(:,2).^2 + E3(:,3).^2;
            
            % Determinant (Sarrus rules)
            detG = a.*e.*i + b.*f.*g + c.*d.*h - c.*e.*g - f.*h.*a - i.*b.*d;
            
            % Inverse by co-factor
            Dx1 = (e.*i - f.*h) ./ detG;
            Dx2 = (c.*h - b.*i) ./ detG;
            Dx3 = (b.*f - c.*e) ./ detG;
            Dy1 = (f.*g - d.*i) ./ detG;
            Dy2 = (a.*i - c.*g) ./ detG;
            Dy3 = (c.*d - a.*f) ./ detG;
            Dz1 = (d.*h - e.*g) ./ detG;
            Dz2 = (b.*g - a.*h) ./ detG;
            Dz3 = (a.*e - b.*d) ./ detG;
            
            % Gradient projection to integration points
            dbas = cell(1,3);
            for n = 1:3
                dbas{n} = zeros(Nelt,Nbas,Ngss);
                DCVx    = Dx1.*E1(:,n) + Dx2.*E2(:,n) + Dx3.*E3(:,n);
                DCVy    = Dy1.*E1(:,n) + Dy2.*E2(:,n) + Dy3.*E3(:,n);
                DCVz    = Dz1.*E1(:,n) + Dz2.*E2(:,n) + Dz3.*E3(:,n);
                for i = 1:Nbas
                    for j = 1:Ngss
                        dbas{n}(:,i,j) = DCVx(:)*dxF(i,j) + DCVy(:)*dyF(i,j) + + DCVz(:)*dzF(i,j);
                    end
                end
            end
        end
        
        % Dof to quadrature matrix
        M = cell(1,3);
        for n = 1:3
            M{n} = sparse(idx(:),jdx(:),dbas{n}(:),Nqud,Ndof);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%% NORMALS * DQM
elseif strcmp(fe.opr,'n*[psi]')
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = '[psi]';
    uqm     = tmp.uqm(domain);
    
    % Normals
    nrm = domain.qudNrm;
    
    % Dot product
    m    = size(nrm,1);
    M{1} = spdiags(nrm(:,1),0,m,m) * uqm;
    M{2} = spdiags(nrm(:,2),0,m,m) * uqm;
    M{3} = spdiags(nrm(:,3),0,m,m) * uqm;


%%% QUADRATURE * DQM
elseif strcmp(fe.opr,'q*[psi]')
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = '[psi]';
    uqm     = tmp.uqm(domain);
    
    % Quadrature
    qud = domain.qud;
    
    % Dot product
    m    = size(qud,1);
    M{1} = spdiags(qud(:,1),0,m,m) * uqm;
    M{2} = spdiags(qud(:,2),0,m,m) * uqm;
    M{3} = spdiags(qud(:,3),0,m,m) * uqm;

    
%%% QUADRATURE dot NORMALS * DQM
elseif strcmp(fe.opr,'qdotn*[psi]')
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = '[psi]';
    uqm     = tmp.uqm(domain);
    
    % Quadrature and normals
    qud = domain.qud;
    nrm = domain.qudNrm;
    
    % Dot product
    m = size(nrm,1);
    M = spdiags(sum(qud.*nrm,2),0,m,m) * uqm;

    
%%%% NORMAL x DQM
elseif strcmp(fe.opr,'nxgrad[psi]')
    % Finite elements
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = 'grad[psi]';
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
      
    
%%% GRADIENT (j)
elseif strcmp(fe.opr(1:end-1),'grad[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite elements
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = 'grad[psi]';
    uqm     = tmp.uqm(domain);
    
    % Gradient (j)
    M  = uqm{j};

        
%%% NORMALS * DQM (j)
elseif strcmp(fe.opr(1:end-1),'n*[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = '[psi]';
    uqm     = tmp.uqm(domain);
    
    % Normal (component j)
    N = domain.qudNrm;
    
    % Dot product (j)
    m = size(N,1);
    M = spdiags(N(:,j),0,m,m) * uqm;
    
    
%%% QUADRATURE * DQM (j)
elseif strcmp(fe.opr(1:end-1),'q*[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = '[psi]';
    uqm     = tmp.uqm(domain);
    
    % Quadrature (component j)
    qud = domain.qud;
    
    % Dot product (j)
    m = size(qud,1);
    M = spdiags(qud(:,j),0,m,m) * uqm;
    
    
%%%% NORMAL x DQM (j)
elseif strcmp(fe.opr(1:end-1),'nxgrad[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite element
    tmp     = fem(fe.msh,fe.typ);
    tmp.opr = 'grad[psi]';
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
    error('femLagrangePn.m : unavailable case')
end
end
