classdef auxiliaryFunctions < handle

    methods (Static, Access = public)
        %% createBoundaryConditions
        function createBoundaryConditions(obj)            
            obj.bcU = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1 2]);
            obj.bcT = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1 2]);
            obj.bcW = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1]);            
        end
        %% createLHS
        function LHS = createLHS(obj)
            % E = obj.young;
            % A = obj.area;
            % f = @(u,v) E.*A.*DDP(SymGrad(v),SymGrad(u));
            
            A = obj.A_coupling;
            f = @(u,v) A.*DDP(SymGrad(v),SymGrad(u));

            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = reduceMatrix(obj,Ku,obj.bcU,obj.bcU);

            % E = obj.young;
            % I = obj.inertia;
            % f = @(u,v) E.*I.*DDP(SymGrad(v),SymGrad(u));
            
            D = obj.D_bending;
            f = @(u,v) D.*DDP(SymGrad(v),SymGrad(u));
            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);


            % A = obj.area;
            % G = obj.shear;
            % f = @(u,v) G.*A.*DP(v,u);

            H = obj.H_shear; 
            f = @(u,v) H.*DP(v,u);
            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);


            A = obj.area;
            G = obj.shear;
            f = @(u,v) H.*DP(v,Grad(u));
            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);            
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);

        
            A = obj.area;
            G = obj.shear;
            f = @(u,v) H.*DP(Grad(v),Grad(u));
            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);  
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            beta = 1; 

            Zut = zeros(nU,nTheta);
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Zut' (Ktheta+beta*Mtheta) beta*Nthetaw; Zuw' beta*Nthetaw' beta*Kw];
        end

        %% createRHS
        function RHS = createRHS(obj)

            p = ConstantFunction.create([0 0],obj.mesh);
            m = ConstantFunction.create([0 0],obj.mesh);
            q = ConstantFunction.create(0,obj.mesh);   

            fu = @(v) DP(p,v);
            RHSu = IntegrateRHS(fu,obj.uFun,obj.mesh,'Domain',2);
            RHSu = obj.reduceVector(RHSu,obj.bcU);

            ftheta   = @(v) DP(m,v);
            RHStheta = IntegrateRHS(ftheta,obj.thetaFun,obj.mesh,'Domain',2);
            RHStheta = obj.reduceVector(RHStheta,obj.bcT);

            fw = @(v) q.*v;
            RHSw = IntegrateRHS(fw,obj.wFun,obj.mesh,'Domain',2);
            RHSw = obj.reduceVector(RHSw,obj.bcW);

            if obj.bcCase == 1
                auxiliaryFunctions.computeForces(obj);
                RHSq = obj.RHSq;
                RHSq = obj.reduceVector(RHSq,obj.bcW);
                RHSw = RHSw + RHSq;
            end            

            RHS = [RHSu;RHStheta;RHSw];
        end


        %% computeForces
        function computeForces(obj)
            bc  = obj.bcW;
            t   = bc.tractionFun;
            rhs = zeros(obj.wFun.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(obj.wFun);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.bcW;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    % Convertir DOFs locales de w a globales (solo libres)
                    freeDofsW = computeFreeDofs(obj,obj.bcW);
                    globalDofsW = auxiliaryFunctions.localToGlobalDofs(obj,freeDofsW, 'w');

                    R = -obj.lhs(globalDofsW, dirich)*dirichV;
                else
                    R = zeros(sum(obj.wFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            obj.RHSq = rhs;
        end

        %% localToGlobalDofs
        function globalDofs = localToGlobalDofs(obj, localDofs, field)
            % Convierte DOFs locales a DOFs globales en la matriz reducida
            % field puede ser 'u', 'theta' o 'w'

            nU = length(obj.bcU.free_dofs);
            nTheta = length(obj.bcT.free_dofs); 

            switch field
                case 'u'
                    offset = 0;
                case 'theta'
                    offset = nU;
                case 'w'
                    offset = nU + nTheta;
                otherwise
                    error('Field must be u, theta or w');
            end

            % Los localDofs ya son índices dentro del conjunto de DOFs libres
            % Solo necesitamos sumarles el offset
            globalDofs = length(localDofs) + offset;
        end


        %% createGeneralBoundaryConditions
        function bc = createGeneralBoundaryConditions(obj,direct)
            TOL = 1e-12;
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;
            isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;
            isTop   = @(coor)  abs(coor(:,2)-yMax)< TOL;
            isBotom = @(coor)  abs(coor(:,2)-yMin)< TOL;


            % Boundary conditions Case
            obj.bcCase = 1;     % 1 --> Modified bc // 2 --> Original bc (change q) 
            
            % APPLYED FORCE
            applyedForce = 1;

            % 1 --> Single node
            % 2 --> Right edge

            switch obj.bcCase
                case 1

                    sDir{1}.domain    = @(coor) isLeft(coor);
                    sDir{1}.direction = direct;
                    sDir{1}.value     = 0;
                    sDir{1}.ndim      = length(direct);


                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;

                    

                    switch applyedForce
                        case 1

                            % Apply point load to a single node on the right boundary
                            % Find right boundary nodes
                            rightNodes = find(isRight(obj.mesh.coord));
                            if isempty(rightNodes)
                                error('No nodes found on the right boundary.');
                            end
                            % Choose one node: here pick the middle one (can change as needed)
                            idx = ceil(numel(rightNodes)/2);
                            singleNode = rightNodes(idx);

                            % Define point load structure for that single node
                            sPL{1}.domain    = @(coor) (1:size(coor,1))'==singleNode;
                            sPL{1}.direction = 2;
                            sPL{1}.value     = 1;

                        case 2

                            % Load on Right edge
                            sPL{1}.domain    = @(coor) isRight(coor);
                            sPL{1}.direction = 2;
                            sPL{1}.value     = 1;
                    end

                    pointloadFun = [];
                    for i = 1:numel(sPL)
                        pl = TractionLoad2(obj.mesh, sPL{i}, 'DIRAC');
                        pointloadFun = [pointloadFun, pl];

                    end
                    s.pointloadFun = pointloadFun;


                    s.periodicFun  = [];
                    s.mesh = obj.mesh;
                    bc = BoundaryConditions(s);

                case 2

                    sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isBotom(coor) | isTop(coor);
                    sDir{1}.direction = direct;
                    sDir{1}.value     = 0;
                    sDir{1}.ndim      = length(direct);

                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;

                    s.pointloadFun = [];
                    s.periodicFun  = [];
                    s.mesh = obj.mesh;
                    bc = BoundaryConditions(s);

            end

        end

        %% createStrainStressFunctions_InternalForces
        function [strainFun, stressFun, Nfun, Mfun, Qfun] = createStrainStressFunctions_InternalForces(obj, z)

            nNodes = size(obj.mesh.coord, 1);

            % Compute layer interfaces from middle plane
            zInterfaces = zeros(length(obj.area.constant)+1, 1);
            zInterfaces(1) = -sum(obj.area.constant)/2;
            for i = 2:length(zInterfaces)
                zInterfaces(i) = zInterfaces(i-1) + obj.area.constant(i-1);
            end

            % Compute SymGrad of u and theta outside the loop (independent of z)
            strainU     = SymGrad(obj.uFun);
            strainTheta = SymGrad(obj.thetaFun);

            strainUvalues     = strainU.evaluate(obj.mesh);
            strainThetaValues = strainTheta.evaluate(obj.mesh);

            nGauss = size(strainUvalues, 3);
            nElem  = size(obj.mesh.connec, 1);

            % Compute gradW at Gauss points and project to nodes
            gradW       = Grad(obj.wFun);
            gradWvalues = gradW.evaluate(obj.mesh);
            dw_dx_gauss = squeeze(gradWvalues(1,:,:));  % (nGauss x nElem)
            dw_dy_gauss = squeeze(gradWvalues(2,:,:));  % (nGauss x nElem)

            % Project gradW to nodes
            dw_dx_nodal = zeros(nNodes, 1);
            dw_dy_nodal = zeros(nNodes, 1);
            node_count  = zeros(nNodes, 1);

            % Project epsilonU to nodes
            epsilonU_nodal     = zeros(nNodes, 3);  % [exx, eyy, exy]
            epsilonTheta_nodal = zeros(nNodes, 3);

            for e = 1:nElem
                elemNodes = obj.mesh.connec(e, :);
                for n = 1:length(elemNodes)
                    nodeIdx = elemNodes(n);
                    dw_dx_nodal(nodeIdx) = dw_dx_nodal(nodeIdx) + mean(dw_dx_gauss(e));
                    dw_dy_nodal(nodeIdx) = dw_dy_nodal(nodeIdx) + mean(dw_dy_gauss(e));

                    epsilonU_nodal(nodeIdx, 1) = epsilonU_nodal(nodeIdx, 1) + mean(squeeze(strainUvalues(1,1,:,e)));
                    epsilonU_nodal(nodeIdx, 2) = epsilonU_nodal(nodeIdx, 2) + mean(squeeze(strainUvalues(2,2,:,e)));
                    epsilonU_nodal(nodeIdx, 3) = epsilonU_nodal(nodeIdx, 3) + mean(squeeze(strainUvalues(1,2,:,e)));

                    epsilonTheta_nodal(nodeIdx, 1) = epsilonTheta_nodal(nodeIdx, 1) + mean(squeeze(strainThetaValues(1,1,:,e)));
                    epsilonTheta_nodal(nodeIdx, 2) = epsilonTheta_nodal(nodeIdx, 2) + mean(squeeze(strainThetaValues(2,2,:,e)));
                    epsilonTheta_nodal(nodeIdx, 3) = epsilonTheta_nodal(nodeIdx, 3) + mean(squeeze(strainThetaValues(1,2,:,e)));

                    node_count(nodeIdx)  = node_count(nodeIdx) + 1;
                end
            end

            dw_dx_nodal = dw_dx_nodal ./ node_count;
            dw_dy_nodal = dw_dy_nodal ./ node_count;

            epsilonU_nodal     = (epsilonU_nodal     ./ node_count).';
            epsilonTheta_nodal = (epsilonTheta_nodal ./ node_count).';

            % Transverse shear strains at nodes
            epsilon_xz_nodal = 0.5 * (dw_dx_nodal + obj.thetaFun.fValues(:,1));
            epsilon_yz_nodal = 0.5 * (dw_dy_nodal + obj.thetaFun.fValues(:,2));

            for i = 1:numel(z)
                % Obtain membrane strains at Gauss points: epsilon = epsilonU + z*epsilonTheta
                epsilon_xx_gauss = squeeze(strainUvalues(1,1,:,:) + z{i} * strainThetaValues(1,1,:,:));
                epsilon_yy_gauss = squeeze(strainUvalues(2,2,:,:) + z{i} * strainThetaValues(2,2,:,:));
                epsilon_xy_gauss = squeeze(strainUvalues(1,2,:,:) + z{i} * strainThetaValues(1,2,:,:));

                % Project membrane strains to nodes
                strain_nodal = zeros(nNodes, 6);
                node_count   = zeros(nNodes, 1);

                for e = 1:nElem
                    elemNodes = obj.mesh.connec(e, :);
                    for n = 1:length(elemNodes)
                        nodeIdx = elemNodes(n);
                        strain_nodal(nodeIdx, 1) = strain_nodal(nodeIdx, 1) + mean(epsilon_xx_gauss(e));
                        strain_nodal(nodeIdx, 2) = strain_nodal(nodeIdx, 2) + mean(epsilon_yy_gauss(e));
                        strain_nodal(nodeIdx, 3) = 0;  % epsilon_zz = 0
                        strain_nodal(nodeIdx, 4) = strain_nodal(nodeIdx, 4) + mean(epsilon_xy_gauss(e));
                        node_count(nodeIdx) = node_count(nodeIdx) + 1;
                    end
                end

                strain_nodal(:,1:4) = strain_nodal(:,1:4) ./ node_count;

                % Add transverse shear strains (already at nodes)
                strain_nodal(:, 5) = epsilon_xz_nodal;
                strain_nodal(:, 6) = epsilon_yz_nodal;

                % Create LagrangianFunction for strains [exx, eyy, ezz, exy, exz, eyz]
                strainFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                strainFun{i}.setFValues(strain_nodal);

                % Obtain constitutive matrix for the layer at z{i}
                C = auxiliaryFunctions.createConstitutiveMatrix(obj, z{i}, zInterfaces);

                % Obtain stress: sigma = C * epsilon at each node
                stress_nodal = (C * strain_nodal.').';  % (nNodes x 6)

                % Create LagrangianFunction for stresses [sxx, syy, szz, sxy, sxz, syz]
                stressFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                stressFun{i}.setFValues(stress_nodal);
            end

            A_Matrix = zeros(3);
            B_Matrix = zeros(3);
            D_Matrix = zeros(3);
            H_Matrix = zeros(2);

            % Compute layer mid-surface positions using numeric arrays instead of cell arrays
            nLayers = numel(obj.zLayer);
            h = zeros(1, nLayers);
            for i = 1:nLayers
                h(i) = obj.zLayer{i}; %- sum(obj.area.constant)/2
            end

            % Layer medium-plane (z) as a numeric vector
            z = 0.5 * (h(2:end) + h(1:end-1));


            for kappa = 1:length(obj.young.constant)
                ConstitutiveMatrix = auxiliaryFunctions.createConstitutiveMatrix(obj,z(kappa),zInterfaces);
                A_Matrix = A_Matrix + ConstitutiveMatrix([1,2,4],[1,2,4])*obj.area.constant(kappa);
                B_Matrix = B_Matrix + 0.5*ConstitutiveMatrix([1,2,4],[1,2,4])*(obj.zLayer{kappa+1}^2-obj.zLayer{kappa}^2);
                D_Matrix = D_Matrix + 1/3*ConstitutiveMatrix([1,2,4],[1,2,4])*(obj.zLayer{kappa+1}^3-obj.zLayer{kappa}^3);
                H_Matrix = H_Matrix + ConstitutiveMatrix(5:6,5:6).*obj.area.constant(kappa);

            end

            Nvalues = A_Matrix*epsilonU_nodal + B_Matrix*epsilonTheta_nodal;
            Mvalues = B_Matrix*epsilonU_nodal + D_Matrix*epsilonTheta_nodal;
            Qvalues = H_Matrix*[dw_dx_nodal.' ; dw_dy_nodal.'];

            Nfun = LagrangianFunction.create(obj.mesh,length(Nvalues(:,1)),'P1');
            Mfun = LagrangianFunction.create(obj.mesh,length(Mvalues(:,1)),'P1');
            Qfun = LagrangianFunction.create(obj.mesh,length(Qvalues(:,1)),'P1');
            Nfun.setFValues(Nvalues.');
            Mfun.setFValues(Mvalues.');
            Qfun.setFValues(Qvalues.');

        end


        %% createConstitutiveMatrix
        function C = createConstitutiveMatrix(obj,z,z_interfaces)
            switch obj.materialProperties
                case 'ISOTROPIC'
                    kLayer = find(z >= z_interfaces(1:(end-1)) & z <= z_interfaces(2:end), 1);
                    nu = obj.poisson.constant(kLayer);
                    lambda = obj.young.constant(kLayer) / (1+nu) / (1-2*nu);
                    C = lambda* [1-nu, nu, nu, 0, 0, 0;
                        nu, 1-nu, nu, 0, 0, 0;
                        nu, nu, 1-nu, 0, 0, 0;
                        0, 0, 0, 0.5*(1-2*nu), 0, 0;
                        0, 0, 0, 0, 0.5*(1-2*nu), 0;
                        0, 0, 0, 0, 0, 0.5*(1-2*nu)];


                case 'ANISOTROPIC'
                    switch obj.materialLayers
                        case 'SINGLE'

                        case 'MULTI'
                    end
            end
        end



    end



end


