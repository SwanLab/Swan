classdef PhaseFieldComputer < handle

    properties (Access = public)

    end

    properties (Access = private)
        mesh
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
    end

    methods (Access = public)

        function obj = PhaseFieldComputer()
            close all %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.init()
            obj.computeMesh();
            obj.createBoundaryConditions();
            obj.createPhaseField();
            obj.createMaterialInterpolation();
            obj.createDissipationInterpolation();

            for i = 1:5
            obj.phaseField.plot;
            title('Phase Field')
            obj.computeFEM();

            obj.createEnergyMassMatrix();
            obj.createDissipationMassMatrix();
            obj.createStiffnessMatrix();
            obj.createEnergyForceVector();
            obj.createDissipationForceVector();
            obj.createForceDerivativeVector();
            obj.solvePhaseFieldEquation();

            obj.phaseField.fValues = obj.phaseField.fValues + obj.deltaPhi;
            end

            obj.phaseField.plot;
            title('Phase Field')
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
            x1 = linspace(0,2,20);
            x2 = linspace(1,2,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;
        end

        function createBoundaryConditions(obj)
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
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
            bc.pointload(:,3) = -0.05;
            obj.boundaryConditions = bc;
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
            %sAF.fHandle = @(x) x(1,:,:)/xmax;
            sAF.fHandle = @(x) x(1,:,:)-x(1,:,:);
            sAF.ndimf   = 1;
            sAF.mesh    = obj.mesh;
            xFun = AnalyticalFunction(sAF);

            phi = xFun.project('P1');
            %phi.plot();
            %title('Initial phi')
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

        %% ELASTIC EQUATION 
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
            obj.fem = FEM.create(s);
            obj.fem.solve();

            %obj.fem.uFun.plot()
            %obj.fem.stressFun.plot()
            %obj.fem.strainFun.plot()

            obj.fem.uFun.fValues(:,end+1) = 0;
            obj.fem.uFun.ndimf = 3;
            obj.fem.print('Example','Paraview')
        end

        %% PHASE-FIELD EQUATION (LHS)
        % Internal energy mass matrix
        function DDenergy = createFGaussDDEnergyFunction(obj)
            s.fValues = obj.computeSecondEnergyDerivativeField();
            s.quadrature = obj.fem.strainFun.quadrature;
            s.mesh = obj.mesh;
            DDenergy = FGaussDiscontinuousFunction(s);
            %DDenergy.plot();
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
            
            nGauss = length(e.fValues);
            DDenergyVal = zeros(1,nGauss);
            for iGauss=1:nGauss
                DDenergyVal(iGauss) = e.fValues(:,:,iGauss)'*(DDmat.C(:,:,iGauss)*e.fValues(:,:,iGauss));
            end
        end

        function createEnergyMassMatrix(obj)
            DDenergyFun =  obj.createFGaussDDEnergyFunction(); 
            s.trial = P1Function.create(obj.mesh,1);
            s.test = P1Function.create(obj.mesh,1);
            s.function = DDenergyFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = DDenergyFun.quadrature.order;
            LHS = LHSintegrator.create(s);
            obj.Mi = LHS.compute(); 
        end
        
        % Dissipation mass matrix
        function DDalpha = createFGaussDDDissipationFunction(obj)
            phi0 = obj.phaseField.project('P0');
            s.fValues = reshape(obj.dissipationInterpolation.computeDDAlphaProp(squeeze(phi0.fValues)),1,1,[]);
            s.quadrature = obj.fem.strainFun.quadrature;
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
            s.quadratureOrder = DDalphaFun.quadrature.order;
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
            s.quadrature = obj.fem.strainFun.quadrature;
            s.mesh = obj.mesh;
            Denergy = FGaussDiscontinuousFunction(s);
            %Denergy.plot();
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
            
            nGauss = length(e.fValues);
            DenergyVal = zeros(1,nGauss);
            for iGauss=1:nGauss
                DenergyVal(iGauss) = e.fValues(:,:,iGauss)'*(Dmat.C(:,:,iGauss)*e.fValues(:,:,iGauss));
            end
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
            s.quadrature = obj.fem.strainFun.quadrature;
            s.mesh = obj.mesh;
            Dalpha = FGaussDiscontinuousFunction(s);
        end   

        function createDissipationForceVector(obj)
            DalphaFun =  obj.createFGaussDDissipationFunction(); 
            test = P1Function.create(obj.mesh,1);
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadratureOrder = DalphaFun.quadrature.order;
            RHS = RHSintegrator.create(s);
            obj.Fd = RHS.compute(DalphaFun, test);     
        end
       
        % Force derivative vector
        function createForceDerivativeVector(obj)
            quad = obj.fem.strainFun.quadrature;
            PhiGradient = obj.phaseField.computeGradient(quad);
            test = P1Function.create(obj.mesh,2);
            s.quadratureOrder = quad.order;
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            RHS = RHSintegrator.create(s);
            forceVector = RHS.compute(PhiGradient, test); 
            obj.DF = forceVector.fValues;
        end

        % Matrix euqation
        function solvePhaseFieldEquation(obj)
            LHS = obj.Mi + obj.Md + obj.K;
            RHS = -(obj.Fi + obj.Fd + obj.DF);
            obj.deltaPhi = LHS\RHS;
        end

    end

end
