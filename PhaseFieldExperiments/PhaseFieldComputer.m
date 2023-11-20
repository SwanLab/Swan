classdef PhaseFieldComputer < handle

    properties (Constant, Access = public)
        tolErrU = 1e-6;
        tolErrPhi = 1e-6;
        Gc = 0.02;
        fc = 1;
        Force = 1;
    end

    properties (Access = private)
        mesh
        quadrature
        boundaryConditions
        materialInterpolation
        material
        dissipationInterpolation
        phaseField
        deltaPhi
        fem
    end

    properties (Access = private)
        Mi
        Md
        K
        Fi
        Fd
        DF
        Constant
        l0
    end

    methods (Access = public)

        function obj = PhaseFieldComputer()
            close all %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.init()
            obj.computeMesh();
            obj.createQuadrature();
            obj.createPhaseField();
            obj.createMaterialInterpolation();
            obj.createDissipationInterpolation();

            %obj.l0 = (27/256)*(1*obj.Gc/obj.fc^2);
            obj.Constant = obj.Gc/(4*0.5);

            obj.l0 = 0.001;

            niter = 100;
            Energy = zeros(4,niter);
            ForceDisplacement = zeros(2,niter);
            Iterations = zeros(2,niter);
            F = zeros(3,11);
            for i = 1:niter
                obj.createBoundaryConditions(i,niter);
                errorU = 1;
                Uold = zeros(obj.mesh.nnodes,obj.mesh.ndim);
                numIterU = 1;
                numIterP = 1;
                while (errorU > obj.tolErrU) && (numIterU < 100)
                    obj.computeFEM();
                    obj.fem.uFun.plot()
                    errorPhi = 1;
                    while (errorPhi > obj.tolErrPhi) && (numIterP < 100)
                        obj.createEnergyMassMatrix();
                        obj.createDissipationMassMatrix();
                        obj.createStiffnessMatrix();
                        obj.createEnergyForceVector();
                        obj.createDissipationForceVector();
                        obj.createForceDerivativeVector();
                        obj.solvePhaseFieldEquation();
                        %obj.computeTotalEnergy()

                        obj.phaseField.fValues = obj.phaseField.fValues + obj.deltaPhi;
                        res = obj.computeResidual();
                        errorPhi = norm(res)
                        numIterP = numIterP + 1;
                        %obj.phaseField.plot;  
                    end
                    errorU = norm(obj.fem.uFun.fValues - Uold)
                    Uold = obj.fem.uFun.fValues;
                    numIterU = numIterU + 1;
                end

                F(1,i) = mean(obj.Fi);
                F(2,i) = mean(obj.Fd);
                F(3,i) = mean(obj.DF);

                Energy(1,i) = obj.computeTotalExternalWork();
                Energy(2,i) = obj.computeTotalEnergy();
                Energy(3,i) = obj.computeTotalDissipationLocal();
                Energy(4,i) = obj.computeTotalDissipationNonLocal();
                Iterations(1,i) = numIterU;
                Iterations(2,i) = numIterP;
                ForceDisplacement(1,i) = obj.boundaryConditions.pointload(1,3);
                ForceDisplacement(2,i) = max(abs(obj.fem.uFun.fValues(:,2)));

                %obj.fem.uFun.plot;

                obj.phaseField.plot;                
                colorbar
                clim([0 1])
                title('Final phase field')
                fprintf('%%%%%%%%%% NEXT ITERATION %%%%%%%%%%%')
                %obj.fem.print(['Example',num2str(i)],'GiD')
                %obj.phaseField.print(['Example',num2str(i),'Phi'],'GiD')
            end
            figure
            plot(F(1,:))
            hold on
            plot(F(2,:))
            hold on
            plot(F(3,:))
            legend('Fi','Fd','DF')


            figure
            plot(Energy(1,:))
            hold on
            plot(Energy(2,:))
            hold on
            plot(Energy(3,:))
            hold on
            plot(Energy(4,:))
            legend('External Work','Internal Energy','Local surface energy','Non-local surface energy')

            figure
            plot(ForceDisplacement(2,:),ForceDisplacement(1,:))
            figure
            plot(ForceDisplacement(2,:))

            figure
            plot(Iterations(1,:))
            hold on
            plot(Iterations(2,:))
            legend('U','phi')
        end

    end




    methods (Access = private)

        function init(obj)
        end

        function computeMesh(obj)
            %file = 'test2d_triangle';
            %a.fileName = file;
            % s = FemDataContainer(a);

            % Generate coordinates
            x1 = linspace(0,1,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectiviti            % x1 = linspace(0,1,20);
            x2 = linspace(0,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh(sBg);
            bdMesh  = bgMesh.createBoundaryMesh();

            sLS.type       = 'circleInclusion';
            sLS.mesh       = bgMesh;
            sLS.ndim       = 2;
            sLS.fracRadius = 0.4;
            sLS.coord      = bgMesh.coord;
            ls = LevelSetCreator.create(sLS);
            levelSet = ls.getValue();

            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bdMesh;
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);

            obj.mesh = uMesh.createInnerMesh();
            obj.mesh = obj.mesh.computeCanonicalMesh();
            %obj.mesh = bgMesh;

            figure
            uMesh.plot
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;

            x1 = linspace(0,1,20);
            x2 = linspace(0,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            sBg.coord  = V(:,1:2);
            sBg.connec = F;
            bgMesh = Mesh(sBg);
            bdMesh  = bgMesh.createBoundaryMesh();

            sLS.type       = 'circleInclusion';
            sLS.mesh       = bgMesh;
            sLS.ndim       = 2;
            sLS.fracRadius = 0.4;
            sLS.coord      = bgMesh.coord;
            ls = LevelSetCreator.create(sLS);
            levelSet = ls.getValue();

            sUm.backgroundMesh = bgMesh;
            sUm.boundaryMesh   = bdMesh;
            uMesh = UnfittedMesh(sUm);
            uMesh.compute(levelSet);

            %obj.mesh = uMesh.createInnerMesh();
            obj.mesh = uMesh.createInnerMeshGoodConditioning()
            %obj.mesh = obj.mesh.computeCanonicalMesh();
            obj.mesh.plot;
            %obj.mesh = bgMesh;

            figure
            uMesh.plot

            % sM.coord = [-1,-1;
            %             1,-1;
            %             1,1;
            %             -1,1];
            % sM.connec = [1 2 3 4];
            % m = Mesh(sM);
            % m.plot();
            % obj.mesh = m;

        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createPhaseField(obj)
            % sLS.type       = 'rectangleInclusion';
            % sLS.mesh       = obj.mesh;
            % sLS.ndim       = 2;
            % sLS.widthH = 1;
            % sLS.widthV = 2;
            % sLS.coord      = obj.mesh.coord;
            % ls = LevelSetCreator.create(sLS);
            % dmg = ls.getValue();
            % obj.phi = dmg;
            % scatter3(obj.mesh.coord(:,1),obj.mesh.coord(:,2),obj.phi);

            xmax = max(obj.mesh.coord(:,1));
            sAF.fHandle = @(x) x(1,:,:)-x(1,:,:);
            %sAF.fHandle = @(x) x(1,:,:)/xmax;
            %sAF.fHandle = @(x) (x(1,:,:).^2)/xmax;
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);

            phi = xFun.project('P1');
            obj.phaseField = phi;
        end

        function createMaterialInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'PhaseFieldI';
            c.nElem = obj.mesh.nelem;
            c.dim = '2D';
            c.constitutiveProperties.rho_plus = 1;
            c.constitutiveProperties.rho_minus = 0;
            c.constitutiveProperties.E_plus = 1;
            c.constitutiveProperties.E_minus = 1e-3;
            c.constitutiveProperties.nu_plus = 1/3;
            c.constitutiveProperties.nu_minus = 1/3;

            matInt = MaterialInterpolation.create(c);
            obj.materialInterpolation = matInt;
        end

        function createDissipationInterpolation(obj)
            c.typeOfMaterial = 'ISOTROPIC';
            c.interpolation = 'PhaseFieldD';

            disInt = MaterialInterpolation.create(c);
            obj.dissipationInterpolation = disInt;
        end

        function createBoundaryConditions(obj,Fstep,niter)
            % dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            % rightSide  = max(obj.mesh.coord(:,1));
            % isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            % isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            % forceNodes = isInRight & isInMiddleEdge;
            % nodes = 1:obj.mesh.nnodes;
            % 
            % ndim = 2;
            % bc.dirichlet = zeros(ndim*length(nodes(dirichletNodes)),3);
            % for i=1:ndim
            %     bc.dirichlet(i:2:end,1) = nodes(dirichletNodes);
            %     bc.dirichlet(i:2:end,2) = i;
            % end
            % bc.dirichlet(:,3) = 0;
            % 
            % bc.pointload(:,1) = nodes(forceNodes);
            % bc.pointload(:,2) = 2;
            % bc.pointload(:,3) = -obj.Force*(Fstep/niter);
            % obj.boundaryConditions = bc;

            dirichletNodes = abs(obj.mesh.coord(:,2)-(0)) < 1e-12;
            UpSide  = max(obj.mesh.coord(:,2));
            isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            forceNodes = isInUp;
            nodes = 1:obj.mesh.nnodes;

            ndim = 2;
            bc.dirichlet = zeros(ndim*length(nodes(dirichletNodes)),3);
            for i=1:ndim
                bc.dirichlet(i:2:end,1) = nodes(dirichletNodes);
                bc.dirichlet(i:2:end,2) = i;
            end
            bc.dirichlet(:,3) = 0;

            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = 1;%obj.Force*(Fstep/niter);
            obj.boundaryConditions = bc;

               %%%%%% TRACTION BC %%%%%%%%%%%%%
            % DownSide = min(obj.mesh.coord(:,2));
            % isInDown = abs(obj.mesh.coord(:,2)-DownSide) < 1e-12;
            % UpSide  = max(obj.mesh.coord(:,2));
            % isInUp = abs(obj.mesh.coord(:,2)-UpSide)< 1e-12;
            % forceNodes = isInUp;
            % nodes = 1:obj.mesh.nnodes;
            % 
            % ndim = 2;
            % DirichletDown = zeros(ndim*length(nodes(isInDown)),3);
            % DirichletUp = zeros(ndim*length(nodes(isInUp)),3);
            % for i=1:ndim
            %     DirichletDown(i:2:end,1) = nodes(isInDown);
            %     DirichletDown(i:2:end,2) = i;
            % 
            %     DirichletUp(i:2:end,1) = nodes(isInUp);
            %     DirichletUp(i:2:end,2) = i;
            % end
            % DirichletDown(:,3) = 0;
            % DirichletUp(1:2:end,3) = 1;
            % 
            % bc.dirichlet = [DirichletDown; DirichletUp];
            % 
            % bc.pointload(:,1) = nodes(forceNodes);
            % bc.pointload(:,2) = 2;
            % bc.pointload(:,3) = 0;%obj.Force*(Fstep/niter);
            % obj.boundaryConditions = bc

        end

        %% %%%%%%%%%%%%%%%%%%%%%% ELASTIC EQUATION %%%%%%%%%%%%%%%%%%%%%%%% %%
        function mat = createMaterial(obj)
            phi0 = obj.phaseField.project('P0');
            mat  = obj.materialInterpolation.computeMatProp(squeeze(phi0.fValues));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mat.kappa;
            s.mu    = mat.mu;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end

        function computeFEM(obj)
            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.bc = obj.boundaryConditions;
            s.interpolationType = 'LINEAR';
            obj.fem = FEM.create(s);
            obj.fem.solve();

            %obj.fem.uFun.plot()
            %obj.fem.stressFun.plot()
            %obj.fem.strainFun.plot()

            %obj.fem.uFun.fValues(:,end+1) = 0;
            %obj.fem.uFun.ndimf = 3;
            %obj.fem.print('Example','Paraview')
        end

        %% %%%%%%%%%%%%%%%% PHASE-FIELD EQUATION (LHS) %%%%%%%%%%%%%%%%%%%%%%%% %%

        % Internal energy mass matrix
        function DDenergy = createFGaussDDEnergyFunction(obj)
            s.fValues = obj.computeSecondEnergyDerivativeField();
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            DDenergy = FGaussDiscontinuousFunction(s);
            DDenergy.plot();
            title('DDEnergy')
        end

        function DDenergyVal = computeSecondEnergyDerivativeField(obj)
            e = obj.fem.strainFun;

            phi0 = obj.phaseField.project('P0');
            DDmat  = obj.materialInterpolation.computeDDMatProp(squeeze(phi0.fValues));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = DDmat.ddkappa;
            s.mu    = DDmat.ddmu;
            DDmat = Material.create(s);
            DDmat.compute(s);

            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            DDenergyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eStre = e.fValues(iStre,iGauss,:).*(DDmat.C(iStre,jStre,:,iGauss).*e.fValues(jStre,iGauss,:));
                        DDenergyVal(iGauss,:) = DDenergyVal(iGauss,:) + squeeze(eStre)';
                    end
                end
            end
            DDenergyVal = 0.5*DDenergyVal;
        end

        function createEnergyMassMatrix(obj)
            DDenergyFun =  obj.createFGaussDDEnergyFunction(); 
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = DDenergyFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = obj.quadrature.order;
            LHS = LHSintegrator.create(s);
            obj.Mi = LHS.compute(); 
        end
        
        % Dissipation mass matrix
        function DDalpha = createFGaussDDDissipationFunction(obj)
            phi0 = obj.phaseField.project('P0');
            s.fValues = reshape(obj.dissipationInterpolation.computeDDAlphaProp(squeeze(phi0.fValues)),1,1,[]);
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            DDalpha = FGaussDiscontinuousFunction(s);
        end

        function createDissipationMassMatrix(obj)
            DDalphaFun =  obj.createFGaussDDDissipationFunction(); 
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = DDalphaFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = obj.quadrature.order;
            LHS = LHSintegrator.create(s);
            obj.Md = LHS.compute(); 
        end

        % Stiffness matrix
        function createStiffnessMatrix(obj)
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'StiffnessMatrix';
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();  
        end

        %% PHASE-FIELD EQUATION (LHS)
        % Internal energy force vector
        function Denergy = createFGaussDEnergyFunction(obj)
            s.fValues = obj.computeFirstEnergyDerivativeField();
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            Denergy = FGaussDiscontinuousFunction(s);
            Denergy.plot();
            title('DEnergy')
        end

        function DenergyVal = computeFirstEnergyDerivativeField(obj)
            e = obj.fem.strainFun;

            phi0 = obj.phaseField.project('P0');
            Dmat  = obj.materialInterpolation.computeDMatProp(squeeze(phi0.fValues));
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = Dmat.dkappa;
            s.mu    = Dmat.dmu;
            Dmat = Material.create(s);
            Dmat.compute(s);

            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            DenergyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eStre = e.fValues(iStre,iGauss,:).*(Dmat.C(iStre,jStre,:,iGauss).*e.fValues(jStre,iGauss,:));
                        DenergyVal(1,iGauss,:) = DenergyVal(iGauss,:) + squeeze(eStre)';
                    end
                end
            end
            DenergyVal = 0.5*DenergyVal;
        end

        function createEnergyForceVector(obj)
            DenergyFun =  obj.createFGaussDEnergyFunction(); 
            test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            RHS = RHSintegrator.create(s);
            obj.Fi = RHS.compute(DenergyFun,test); 
        end
        
        % Dissipation force vector
        function Dalpha = createFGaussDDissipationFunction(obj)
            phi0 = obj.phaseField.project('P0');
            s.fValues = reshape(obj.dissipationInterpolation.computeDAlphaProp(squeeze(phi0.fValues)),1,1,[]);
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            Dalpha = FGaussDiscontinuousFunction(s);
        end   

        function createDissipationForceVector(obj)
            DalphaFun =  obj.createFGaussDDissipationFunction(); 
            test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadratureOrder = obj.quadrature.order;
            RHS = RHSintegrator.create(s);
            obj.Fd = RHS.compute(DalphaFun, test);     
        end
       
        % Force derivative vector
        function createForceDerivativeVector(obj)
            PhiGradient = obj.phaseField.computeGradient(obj.quadrature);
            test = P1Function.create(obj.mesh,1);
            s.quadratureOrder = obj.quadrature.order;
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            RHS = RHSintegrator.create(s);
            forceVector = RHS.compute(PhiGradient, test); 
            obj.DF = forceVector.fValues;
        end

        % Matrix euqation
        function solvePhaseFieldEquation(obj)
            LHS = obj.Mi + (obj.Constant/obj.l0)*obj.Md + (obj.Constant*obj.l0)*obj.K;
            RHS = obj.computeResidual();
            obj.deltaPhi = LHS\RHS;
        end

        function res = computeResidual(obj)
            res = -(obj.Fi + (obj.Constant/obj.l0)*obj.Fd + (obj.Constant*obj.l0)*obj.DF);
        end

        function totVal = computeTotalEnergy(obj)
            e = obj.fem.strainFun;
            mat = obj.material;

            nstre = size(e.fValues,1);
            nGauss = size(e.fValues,2);
            nelem = size(e.fValues,3);
            energyVal = zeros(1,nGauss,nelem);
            for iStre = 1:nstre
                for jStre=1:nstre
                    for iGauss=1:nGauss
                        eStre = e.fValues(iStre,iGauss,:).*(mat.C(iStre,jStre,:,iGauss).*e.fValues(jStre,iGauss,:));
                        energyVal(iGauss,:) = energyVal(iGauss,:) + squeeze(eStre)';
                    end
                end
            end
            energyVal = 0.5*energyVal;
            s.fValues = reshape(energyVal,[1,1,nelem]);
            s.quadrature = obj.quadrature;
            s.mesh = obj.mesh;
            energy = FGaussDiscontinuousFunction(s);
            energy.plot()
            title('Energy')
            
            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'InternalEnergy';
            int = Integrator.create(q);
            totVal = int.compute(e,mat.C);
        end

        function totVal = computeTotalDissipationLocal(obj)
            phi0 = obj.phaseField.project('P0');
            aValues = obj.dissipationInterpolation.computeAlphaProp(squeeze(phi0.fValues));
            s.mesh = obj.mesh;
            s.fValues = aValues;
            alpha = P0Function(s);

            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = (obj.Constant/obj.l0)*int.compute(alpha);
        end

        function totVal = computeTotalDissipationNonLocal(obj)
            PhiGradient = obj.phaseField.computeGradient(obj.quadrature);
            GradGrad = sum(PhiGradient.fValues.^2);
            s.mesh = obj.mesh;
            s.fValues = GradGrad;
            GradGradFun = P0Function(s);
            
            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'Function';
            int = Integrator.create(q);
            totVal = 0.5*(obj.Constant*obj.l0)*int.compute(GradGradFun);
        end

        function totVal = computeTotalExternalWork(obj)
            u = obj.fem.uFun;
            forceValues = zeros(size(u.fValues));
            pLoad = obj.boundaryConditions.pointload;
            idx = sub2ind(size(forceValues),pLoad(:,1),pLoad(:,2));
            forceValues(idx) = pLoad(:,3);

            s.mesh = obj.mesh;
            s.fValues = forceValues;
            f = P1Function(s);
            
            q.mesh = obj.mesh;
            q.quadType = 'CONSTANT';
            q.type = 'ScalarProduct';
            int = Integrator.create(q);
            totVal = int.compute(u,f);
        end

    end

end
