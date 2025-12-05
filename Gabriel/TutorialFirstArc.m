classdef TutorialFirstArc < handle

    properties (Access = private)
        mesh
        filterCompliance
        filterPerimeter
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        perimeter
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = TutorialFirstArc()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterCompliance();
            % obj.createFilterPerimeter();
            %obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createBaseDomain();
            obj.createComplianceFromConstiutive();
            % obj.createVolumeFunctional();
            % obj.createComplianceConstraint();
            obj.createVolumeConstraint();
            % obj.createPerimeter();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            % x1      = linspace(0,2,100);
            % x2      = linspace(0,1,50);
            % [xv,yv] = meshgrid(x1,x2);
            % [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            % s.coord  = V(:,1:2);
            % s.connec = F;
            % obj.mesh = Mesh.create(s);
            obj.mesh = TriangleMesh(2,1,100,50);
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            rho        = DesignVariable.create(sD);

            obj.designVariable = rho;
        end

        function createFilterCompliance(obj)
            s.filterType         = 'LUMP';
            s.mesh               = obj.mesh;
            s.trial              = LagrangianFunction.create(obj.mesh,1,'P1');
            f                    = Filter.create(s);
            obj.filterCompliance = f;
        end

        % function createFilterPerimeter(obj)
        %     s.filterType        = 'PDE';
        %     s.boundaryType      = 'Robin';
        %     s.mesh              = obj.mesh;
        %     s.trial             = LagrangianFunction.create(obj.mesh,1,'P1');
        %     f                   = Filter.create(s);
        %     obj.filterPerimeter = f;
        % end

        % function createMaterialInterpolator(obj)
        %     E0   = 1e-3;
        %     nu0  = 1/3;
        %     E1   = 1;
        %     nu1  = 1/3;
        %     ndim = 2;
        % 
        %     matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
        %     matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);
        % 
        %     matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
        %     matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);
        % 
        %     s.typeOfMaterial = 'ISOTROPIC';
        %     s.interpolation  = 'SIMPALL';
        %     s.dim            = '2D';
        %     s.matA = matA;
        %     s.matB = matB;
        % 
        %     m = MaterialInterpolator.create(s);
        %     obj.materialInterpolator = m;
        % end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end
        function c = createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterCompliance;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            
        end
        % function c = createVolumeFunctional(obj)
        %     s.mesh = obj.mesh;
        %     s.uMesh = obj.createBaseDomain;
        %     s.test = LagrangianFunction.create(obj.mesh,1,'P1');
        %     c = VolumeFunctional(s);
        % end



        % function createComplianceConstraint(obj)
        %     s.mesh                       = obj.mesh;
        %     s.filter                     = obj.filterCompliance;
        %     s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
        %     s.material                   = obj.createMaterial();
        %     s.complianceTarget           = 3;
        %     c = ComplianceConstraint(s);
        %     obj.compliance = c;
        % end
        % 
        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        % function createPerimeter(obj)
        %     eOverhmin     = 10; % 10
        %     epsilon       = eOverhmin*obj.mesh.computeMeanCellSize();
        %     s.mesh        = obj.mesh;
        %     s.filter      = obj.filterPerimeter;
        %     s.epsilon     = epsilon;
        %     s.value0      = 6; % external Perimeter
        %     P             = PerimeterFunctional(s);
        %     obj.perimeter = P;
        % end

        function createCost(obj)
            s.shapeFunctions{1} = obj.createCompliance();
            % s.shapeFunctions{2} = obj.createVolumeFunctional();
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh,'Domain');
        end

        function createConstraint(obj)
            % s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

       function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 5;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            % s.monitoring     = true;
            % s.cost           = obj.cost;
            % s.constraint     = obj.constraint;
            % s.designVariable = obj.designVariable;
            % s.maxIter        = 400;
            % s.ub              = 1;
            % s.lb              = 0;
            % % s.tolerance      = 1e-8;
            % s.constraintCase = {'EQUALITY'};
            % % s.primalUpdater  = obj.primalUpdater;
            % % s.etaNorm        = 0.02;
            % % s.etaNormMin     = 0.02;
            % % s.gJFlowRatio    = 1;
            % % s.etaMax         = 1;
            % % s.etaMaxMin      = 0.01;
            % opt = OptimizerNullSpace(s);
            % % opt = OptimizerProjectedGradient(s);
            % opt.solveProblem();
            % obj.optimizer = opt;

            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 800;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 2;
            s.gif            = true;
            s.gifName        = 'Tutorial_Homo_ReinforcedHexagon_Arch';
            s.printing       = true;
            s.printName      = 'Tutorial_Homo_ReinforcedHexagon_Arch';
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;

            
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1'); 
            % f = obj.filterCompliance.compute(f{1},1);            
           % s.type                 = 'DensityBased';
           % s.density              = f;
           % s.materialInterpolator = obj.materialInterpolator;
           % s.dim                  = '2D';
           % s.mesh                 = obj.mesh;
           
        
            s.density  =  f;
            s.type     = 'HomogenizedMicrostructure';
            s.mesh     = obj.mesh;
            s.young    = 1.0;
            s.fileName = 'HomogenizationResultsReinforcedHexagon';
            m = MaterialFactory.create(s);

           % m = Material.create(s);  


        end

        function bc = createBoundaryConditions(obj)
            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMin = min(obj.mesh.coord(:,2));
        
            % Base do domínio (y = 0)
            isBottom = @(coor) abs(coor(:,2) - yMin) < 1e-12;
        
            % --- Arch: apoios em 0–0.2 e 1.8–2.0 na base ---
            isDirLeft  = @(coor) isBottom(coor) & ...
                                 coor(:,1) >= xMin        & coor(:,1) <= xMin + 0.2;
            isDirRight = @(coor) isBottom(coor) & ...
                                 coor(:,1) >= xMax - 0.2 & coor(:,1) <= xMax;
        
            % --- Arch: carga no trecho central da base (0.9–1.1) ---
            isForce = @(coor) isBottom(coor) & ...
                              coor(:,1) >= xMin + 0.9 & coor(:,1) <= xMin + 1.1;
        
            % Dirichlet: u = [0,0] nos apoios
            sDir{1}.domain    = @(coor) isDirLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
        
            sDir{2}.domain    = @(coor) isDirRight(coor);
            sDir{2}.direction = [1,2];
            sDir{2}.value     = 0;
        
            % Neumann: f = [0,-1] na parte central
            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;
        
            % Montagem das condições de contorno
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
        
            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;
        
            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end