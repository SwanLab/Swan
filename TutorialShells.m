classdef TutorialShells < handle


    properties (Access = private)
        mesh
        young
        area
        shear
        inertia
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
        lhs,RHSS
        solverType
        type, values, fun
        bcCase
    end

    methods (Access = public)

        function obj = TutorialShells()
            clc; close all; 
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()

            obj.solverType = 'REDUCED';
            
            obj.createBoundaryConditions()  % BOUNDARY CONDITIONS 
            LHS = obj.createLHS();
            obj.lhs = LHS;
            RHS = obj.createRHS();          % BC --> Distributed forces 

            x = LHS\RHS;

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));            
            uF = x(1:nU,1);
            tF = x((nU+1):(nU+nTheta),1);
            wF = x((nU+nTheta+1):(nU+nTheta+nW),1);

            dofFT = obj.computeFreeDofs(obj.bcT);
            dofFU = obj.computeFreeDofs(obj.bcU);
            dofFW = obj.computeFreeDofs(obj.bcW);

            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFU,1) = uF; 
            uT = reshape(uT,obj.uFun.ndimf,[])';
            obj.uFun.setFValues(uT);

            wT = zeros(obj.wFun.nDofs,1);
            wT(dofFW,1) = wF; 
            wT = reshape(wT,[], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);
            
            thetaT = zeros(obj.thetaFun.nDofs,1);    
            thetaT(dofFT,1) = tF; 
            thetaT = reshape(thetaT,obj.thetaFun.ndimf,[])';
            obj.thetaFun.setFValues(thetaT);

            z{1} = 0;
            z{2} = 0.1/2;

            [strainFun, stressFun] = obj.createStrainStressFunctions(z);

            %% PLOT AND PRINT 

            plotMatlab = true;
            printParaview = false;

            if plotMatlab == true 
                plot(obj.uFun)
                plot(obj.wFun)
                plot(obj.thetaFun)

                kappa = 1;
                plot(strainFun{kappa});
                title('Epsilon xx')
                plot(stressFun{kappa});
                
            end

            if printParaview == true
                obj.wFun.print('wfun print','Paraview')
                obj.uFun.print('ufun print','Paraview')
                obj.thetaFun.print('thetafun print','Paraview')
                strainFun{kappa}.print('strain', 'Paraview'); % Kappa defined on plots
                stressFun{kappa}.print('stress', 'Paraview');

            end
        end

    end

    methods (Access = private)

        function createMesh(obj)
          obj.mesh = UnitTriangleMesh(5,5);
        end

        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function createMaterialProperties(obj)
          E = 3;
          obj.young = ConstantFunction.create(E,obj.mesh);
          obj.area = ConstantFunction.create(1,obj.mesh);
          obj.shear = ConstantFunction.create(1,obj.mesh);
          obj.inertia = ConstantFunction.create(1,obj.mesh);
        end

        function LHS = createLHS(obj)
            E = obj.young;
            A = obj.area;
            f = @(u,v) E.*A.*DDP(SymGrad(v),SymGrad(u));
            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = obj.reduceMatrix(Ku,obj.bcU,obj.bcU);

            E = obj.young;
            I = obj.inertia;
            f = @(u,v) E.*I.*DDP(SymGrad(v),SymGrad(u));
            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);


            A = obj.area;
            G = obj.shear;
            f = @(u,v) G.*A.*DP(v,u);
            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);


            A = obj.area;
            G = obj.shear;
            f = @(u,v) G.*A.*DP(v,Grad(u));
            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);            
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);

        
            A = obj.area;
            G = obj.shear;
            f = @(u,v) A.*G.*DP(Grad(v),Grad(u));
            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);  
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            Zut = zeros(nU,nTheta);
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; Zut' (Ktheta+Mtheta) Nthetaw; Zuw' Nthetaw' Kw];
        end

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
                obj.computeForces();
                RHSq = obj.RHSS;
                RHSq = obj.reduceVector(RHSq,obj.bcW);
                RHSw = RHSw + RHSq;
            end            

            RHS = [RHSu;RHStheta;RHSw];
        end

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
                    freeDofsW = obj.computeFreeDofs(obj.bcW);
                    globalDofsW = obj.localToGlobalDofs(freeDofsW, 'w');

                    R = -obj.lhs(globalDofsW, dirich)*dirichV;
                else
                    R = zeros(sum(obj.wFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            obj.RHSS = rhs;
        end

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

        function createBoundaryConditions(obj)            
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);            
        end


        % ============================================================================================================================================================
        % ============================================================================================================================================================

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


            % MODIFICATIONS: Embedded beam on one side ant with punctual load on the other  

            obj.bcCase = 1;  % 1 --> Modified bc // 2 --> Original bc (change q) 

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

                    %% APPLYED FORCE 

                    applyedForce = 1;

                    % 1 --> Single node
                    % 2 --> Right edge 

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
        


        % ============================================================================================================================================================
        % ============================================================================================================================================================


        function fD = computeFreeDofs(obj,bC)
            dofs = 1:bC.dirichletFun.nDofs;
            fD   = setdiff(dofs,bC.dirichlet_dofs); 
        end

        function LHSred = reduceMatrix(obj,LHS,bcV,bcU)
            fdofV = obj.computeFreeDofs(bcV);
            fdofU = obj.computeFreeDofs(bcU);
            LHSred = LHS(fdofV,fdofU);
        end

        function RHSred = reduceVector(obj,RHS,bc)
            fdofV = obj.computeFreeDofs(bc);
            RHSred = RHS(fdofV,1);
        end

        function [strainFun, stressFun] = createStrainStressFunctions(obj, z)

            nNodes = size(obj.mesh.coord, 1);

            nu = 0.3;
            lambda = obj.young.constant / (1+nu) / (1-2*nu);

            C = lambda* [1-nu, nu, nu, 0, 0, 0;
                nu, 1-nu, nu, 0, 0, 0;
                nu, nu, 1-nu, 0, 0, 0;
                0, 0, 0, 0.5*(1-2*nu), 0, 0;
                0, 0, 0, 0, 0.5*(1-2*nu), 0;
                0, 0, 0, 0, 0, 0.5*(1-2*nu)];

            for i = 1:numel(z)
                % Obtain strain
                up = obj.uFun.fValues + z{i} * obj.thetaFun.fValues;
                uField = LagrangianFunction.create(obj.mesh,obj.uFun.ndimf,'P1');
                uField.setFValues(up);
                strainDomainFun = SymGrad(uField);

                coords = obj.mesh.coord;

                strainValues = strainDomainFun.evaluate(coords);
                epsilon_xx = squeeze(strainValues(1,1,:,:));
                epsilon_yy = squeeze(strainValues(2,2,:,:));
                epsilon_xy = squeeze(strainValues(1,2,:,:));

                gradW = Grad(obj.wFun);
                gradWvalues = gradW.evaluate(coords);
                dw_dx = squeeze(gradWvalues(1,:,:));
                dw_dy = squeeze(gradWvalues(2,:,:));

                thetaValues = obj.thetaFun.evaluate(coords);
                theta_x = squeeze(thetaValues(1,:,:));
                theta_y = squeeze(thetaValues(2,:,:));

                epsilon_xz = 0.5 * (dw_dx + theta_x);
                epsilon_yz = 0.5 * (dw_dy + theta_y);
                epsilon_zz = zeros(size(epsilon_xx));

                strain = [epsilon_xx(:).'; epsilon_yy(:).'; epsilon_zz(:).'; epsilon_xy(:).'; epsilon_xz(:).'; epsilon_yz(:).'];

                % Obtain stress
                stress = C * strain;

                % Obtain dimensions
                nGauss = size(strainValues,3);
                nElem = size(obj.mesh.connec, 1);

                % Initialize nodal values for strains (nNodes x 6)
                strain_nodal = zeros(nNodes, 6);
                node_count = zeros(nNodes, 1);

                % Reshape strain components
                strain_xx_gauss = reshape(strain(1,:), nGauss, nElem);
                strain_yy_gauss = reshape(strain(2,:), nGauss, nElem);
                strain_zz_gauss = reshape(strain(3,:), nGauss, nElem);
                strain_xy_gauss = reshape(strain(4,:), nGauss, nElem);
                strain_xz_gauss = reshape(strain(5,:), nGauss, nElem);
                strain_yz_gauss = reshape(strain(6,:), nGauss, nElem);

                % Project Gauss points to nodes
                for e = 1:nElem
                    elemNodes = obj.mesh.connec(e, :);
                    for n = 1:length(elemNodes)
                        nodeIdx = elemNodes(n);
                        strain_nodal(nodeIdx, 1) = strain_nodal(nodeIdx, 1) + mean(strain_xx_gauss(:, e));
                        strain_nodal(nodeIdx, 2) = strain_nodal(nodeIdx, 2) + mean(strain_yy_gauss(:, e));
                        strain_nodal(nodeIdx, 3) = strain_nodal(nodeIdx, 3) + mean(strain_zz_gauss(:, e));
                        strain_nodal(nodeIdx, 4) = strain_nodal(nodeIdx, 4) + mean(strain_xy_gauss(:, e));
                        strain_nodal(nodeIdx, 5) = strain_nodal(nodeIdx, 5) + mean(strain_xz_gauss(:, e));
                        strain_nodal(nodeIdx, 6) = strain_nodal(nodeIdx, 6) + mean(strain_yz_gauss(:, e));
                        node_count(nodeIdx) = node_count(nodeIdx) + 1;
                    end
                end

                % Normalize all components
                strain_nodal = strain_nodal ./ node_count;

                % Create single LagrangianFunction for strains with 6 components
                % Order: [exx, eyy, ezz, exy, exz, eyz]
                strainFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                strainFun{i}.setFValues(strain_nodal);

                % Initialize nodal values for stress (nNodes x 6)
                stress_nodal = zeros(nNodes, 6);
                node_count = zeros(nNodes, 1);

                % Reshape stress components
                stress_xx_gauss = reshape(stress(1,:), nGauss, nElem);
                stress_yy_gauss = reshape(stress(2,:), nGauss, nElem);
                stress_zz_gauss = reshape(stress(3,:), nGauss, nElem);
                stress_xy_gauss = reshape(stress(4,:), nGauss, nElem);
                stress_xz_gauss = reshape(stress(5,:), nGauss, nElem);
                stress_yz_gauss = reshape(stress(6,:), nGauss, nElem);

                % Project Gauss points to nodes
                for e = 1:nElem
                    elemNodes = obj.mesh.connec(e, :);
                    for n = 1:length(elemNodes)
                        nodeIdx = elemNodes(n);
                        stress_nodal(nodeIdx, 1) = stress_nodal(nodeIdx, 1) + mean(stress_xx_gauss(:, e));
                        stress_nodal(nodeIdx, 2) = stress_nodal(nodeIdx, 2) + mean(stress_yy_gauss(:, e));
                        stress_nodal(nodeIdx, 3) = stress_nodal(nodeIdx, 3) + mean(stress_zz_gauss(:, e));
                        stress_nodal(nodeIdx, 4) = stress_nodal(nodeIdx, 4) + mean(stress_xy_gauss(:, e));
                        stress_nodal(nodeIdx, 5) = stress_nodal(nodeIdx, 5) + mean(stress_xz_gauss(:, e));
                        stress_nodal(nodeIdx, 6) = stress_nodal(nodeIdx, 6) + mean(stress_yz_gauss(:, e));
                        node_count(nodeIdx) = node_count(nodeIdx) + 1;
                    end
                end

                % Normalize all components
                stress_nodal = stress_nodal ./ node_count;

                % Create single LagrangianFunction for stresses with 6 components
                % Order: [sxx, syy, szz, sxy, sxz, syz]
                stressFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
                stressFun{i}.setFValues(stress_nodal);
            end
        end

        function a = getStrainStress2(obj, z)

        end



    end

end