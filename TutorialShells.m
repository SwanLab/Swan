classdef TutorialShells < handle


    properties (Access = {?TutorialShells, ?auxiliaryFunctions})  %(Access = private)
        mesh
        young
        area
        shear
        inertia
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
        lhs,RHSq
        solverType
        type, values, fun
        bcCase
        aux
    end

    methods (Access = public)

        function obj = TutorialShells()
            % Añadir BETA = 1 
            clc; close all; 

            obj.aux = auxiliaryFunctions;
            
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()

            obj.solverType = 'REDUCED';
            
            obj.createBoundaryConditions()  % BOUNDARY CONDITIONS 
            LHS = auxiliaryFunctions.createLHS(obj);
            obj.lhs = LHS;
            RHS = auxiliaryFunctions.createRHS(obj);          % BC --> Distributed forces 

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

            z{2} = 0;
            z{1} = 0.1/2;

            [strainFun, stressFun] = auxiliaryFunctions.createStrainStressFunctions(obj,z);


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

    methods (Access = {?TutorialShells, ?auxiliaryFunctions})  % (Access = private)

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


        function createBoundaryConditions(obj)            
            obj.bcU = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1 2]);
            obj.bcT = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1 2]);
            obj.bcW = auxiliaryFunctions.createGeneralBoundaryConditions(obj,[1]);            
        end

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


        % 
        % function [strainFun, stressFun] = createStrainStressFunctions(obj, z)
        % 
        %     nNodes = size(obj.mesh.coord, 1);
        % 
        %     nu = 0.3;
        %     lambda = obj.young.constant / (1+nu) / (1-2*nu);
        % 
        %     C = lambda* [1-nu, nu, nu, 0, 0, 0;
        %         nu, 1-nu, nu, 0, 0, 0;
        %         nu, nu, 1-nu, 0, 0, 0;
        %         0, 0, 0, 0.5*(1-2*nu), 0, 0;
        %         0, 0, 0, 0, 0.5*(1-2*nu), 0;
        %         0, 0, 0, 0, 0, 0.5*(1-2*nu)];
        % 
        %     for i = 1:numel(z)
        %         % Obtain strain
        %         up = obj.uFun.fValues + z{i} * obj.thetaFun.fValues;
        %         uField = LagrangianFunction.create(obj.mesh,obj.uFun.ndimf,'P1');
        %         uField.setFValues(up);
        %         strainDomainFun = SymGrad(uField);
        % 
        %         coords = obj.mesh.coord;
        % 
        %         strainValues = strainDomainFun.evaluate(coords);
        %         epsilon_xx = squeeze(strainValues(1,1,:,:));
        %         epsilon_yy = squeeze(strainValues(2,2,:,:));
        %         epsilon_xy = squeeze(strainValues(1,2,:,:));
        % 
        %         gradW = Grad(obj.wFun);
        %         gradWvalues = gradW.evaluate(coords);
        %         dw_dx = squeeze(gradWvalues(1,:,:));
        %         dw_dy = squeeze(gradWvalues(2,:,:));
        % 
        %         thetaValues = obj.thetaFun.evaluate(coords);
        %         theta_x = squeeze(thetaValues(1,:,:));
        %         theta_y = squeeze(thetaValues(2,:,:));
        %         epsilon_xz = 0.5 * (dw_dx + theta_x);
        %         epsilon_yz = 0.5 * (dw_dy + theta_y);
        % 
        % 
        %         epsilon_zz = zeros(size(epsilon_xx));
        % 
        %         strain = [epsilon_xx(:).'; epsilon_yy(:).'; epsilon_zz(:).'; epsilon_xy(:).'; epsilon_xz(:).'; epsilon_yz(:).'];
        % 
        %         % Obtain dimensions
        %         nGauss = size(strainValues,3);
        %         nElem = size(obj.mesh.connec, 1);
        % 
        %         % Initialize nodal values for strains (nNodes x 6)
        %         strain_nodal = zeros(nNodes, 6);
        %         node_count = zeros(nNodes, 1);
        % 
        %         % Reshape strain components
        %         strain_xx_gauss = reshape(strain(1,:), nGauss, nElem);
        %         strain_yy_gauss = reshape(strain(2,:), nGauss, nElem);
        %         strain_zz_gauss = reshape(strain(3,:), nGauss, nElem);
        %         strain_xy_gauss = reshape(strain(4,:), nGauss, nElem);
        %         strain_xz_gauss = reshape(strain(5,:), nGauss, nElem);
        %         strain_yz_gauss = reshape(strain(6,:), nGauss, nElem);
        % 
        % 
        %         %% second methodology
        %         % strainU = SymGrad(obj.uFun);
        %         % strainTheta = SymGrad(obj.thetaFun);
        %         % strainW = SymGrad(obj.wFun);
        %         %
        %         % strainUvalues = strainU.evaluate(obj.mesh.coord);
        %         % strainThetavalues = strainTheta.evaluate(obj.mesh.coord);
        %         % strainWvalues = strainW.evaluate(obj.mesh.coord);
        %         %
        %         % epsilon_xx = squeeze(strainUvalues(1,1,:,:) + z{i}*strainThetavalues(1,1,:,:));
        %         % epsilon_yy = squeeze(strainUvalues(2,2,:,:) + z{i}*strainThetavalues(2,2,:,:));
        %         % epsilon_xy = squeeze(strainUvalues(1,2,:,:) + z{i}*strainThetavalues(1,2,:,:));
        %         %
        %         % gradW = Grad(obj.wFun);
        %         % gradWvalues = gradW.evaluate(obj.mesh.coord);
        %         % dw_dx = squeeze(gradWvalues(1,:,:));
        %         % dw_dy = squeeze(gradWvalues(2,:,:));
        %         %
        %         % thetaValues = obj.thetaFun.evaluate(obj.mesh.coord);
        %         % theta_x = squeeze(thetaValues(1,:,:));
        %         % theta_y = squeeze(thetaValues(2,:,:));
        %         %
        %         % epsilon_xz = 0.5 * (dw_dx + theta_x);
        %         % epsilon_yz = 0.5 * (dw_dy + theta_y);
        %         % epsilon_zz = zeros(size(epsilon_xx));
        %         %
        %         % strain2 = [epsilon_xx(:).'; epsilon_yy(:).'; epsilon_zz(:).'; epsilon_xy(:).'; epsilon_xz(:).'; epsilon_yz(:).'];
        %         %
        %         % diff = abs(strain-strain2);
        %         % nonZeroValues = diff(diff ~= 0)
        %         %
        %         % strain = strain2;
        % 
        % 
        %         %%
        % 
        %         % Project Gauss points to nodes
        %         for e = 1:nElem
        %             elemNodes = obj.mesh.connec(e, :);
        %             for n = 1:length(elemNodes)
        %                 nodeIdx = elemNodes(n);
        %                 strain_nodal(nodeIdx, 1) = strain_nodal(nodeIdx, 1) + mean(strain_xx_gauss(:, e));
        %                 strain_nodal(nodeIdx, 2) = strain_nodal(nodeIdx, 2) + mean(strain_yy_gauss(:, e));
        %                 strain_nodal(nodeIdx, 3) = strain_nodal(nodeIdx, 3) + mean(strain_zz_gauss(:, e));
        %                 strain_nodal(nodeIdx, 4) = strain_nodal(nodeIdx, 4) + mean(strain_xy_gauss(:, e));
        %                 strain_nodal(nodeIdx, 5) = strain_nodal(nodeIdx, 5) + mean(strain_xz_gauss(:, e));
        %                 strain_nodal(nodeIdx, 6) = strain_nodal(nodeIdx, 6) + mean(strain_yz_gauss(:, e));
        %                 node_count(nodeIdx) = node_count(nodeIdx) + 1;
        %             end
        %         end
        % 
        %         % Normalize all components
        %         strain_nodal = strain_nodal ./ node_count;
        % 
        %         % Create single LagrangianFunction for strains with 6 components
        %         % Order: [exx, eyy, ezz, exy, exz, eyz]
        %         strainFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
        %         strainFun{i}.setFValues(strain_nodal);
        % 
        %         % Obtain stress
        %         stress = C * strain;
        % 
        % 
        %         % Initialize nodal values for stress (nNodes x 6)
        %         stress_nodal = zeros(nNodes, 6);
        %         node_count = zeros(nNodes, 1);
        % 
        %         % Reshape stress components
        %         stress_xx_gauss = reshape(stress(1,:), nGauss, nElem);
        %         stress_yy_gauss = reshape(stress(2,:), nGauss, nElem);
        %         stress_zz_gauss = reshape(stress(3,:), nGauss, nElem);
        %         stress_xy_gauss = reshape(stress(4,:), nGauss, nElem);
        %         stress_xz_gauss = reshape(stress(5,:), nGauss, nElem);
        %         stress_yz_gauss = reshape(stress(6,:), nGauss, nElem);
        % 
        %         % Project Gauss points to nodes
        %         for e = 1:nElem
        %             elemNodes = obj.mesh.connec(e, :);
        %             for n = 1:length(elemNodes)
        %                 nodeIdx = elemNodes(n);
        %                 stress_nodal(nodeIdx, 1) = stress_nodal(nodeIdx, 1) + mean(stress_xx_gauss(:, e));
        %                 stress_nodal(nodeIdx, 2) = stress_nodal(nodeIdx, 2) + mean(stress_yy_gauss(:, e));
        %                 stress_nodal(nodeIdx, 3) = stress_nodal(nodeIdx, 3) + mean(stress_zz_gauss(:, e));
        %                 stress_nodal(nodeIdx, 4) = stress_nodal(nodeIdx, 4) + mean(stress_xy_gauss(:, e));
        %                 stress_nodal(nodeIdx, 5) = stress_nodal(nodeIdx, 5) + mean(stress_xz_gauss(:, e));
        %                 stress_nodal(nodeIdx, 6) = stress_nodal(nodeIdx, 6) + mean(stress_yz_gauss(:, e));
        %                 node_count(nodeIdx) = node_count(nodeIdx) + 1;
        %             end
        %         end
        % 
        %         % Normalize all components
        %         stress_nodal = stress_nodal ./ node_count;
        % 
        %         % Create single LagrangianFunction for stresses with 6 components
        %         % Order: [sxx, syy, szz, sxy, sxz, syz]
        %         stressFun{i} = LagrangianFunction.create(obj.mesh, 6, 'P1');
        %         stressFun{i}.setFValues(stress_nodal);
        %     end
        % end

       


        



    end

end