classdef TopOptTestTutorialBoundFormulation < handle

    properties (Access = private)
        mesh
        filterRhoE
        filterRhoI
        filterRhoD
        filterGradientRhoE
        filterGradientRhoI
        filterGradientRhoD
        designVariable
        materialInterpolator
        physicalProblem
        linearBoundFunction
        complianceRhoE
        complianceRhoI
        complianceRhoD
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestTutorialBoundFormulation()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterRhoE();
            obj.createFilterRhoI();
            obj.createFilterRhoD();
            obj.createFilterAdjointRhoE();
            obj.createFilterAdointRhoI();
            obj.createFilterAdjointRhoD();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createLinearBoundFunction();
            obj.createComplianceBoundConstraintRhoE();
            obj.createComplianceBoundConstraintRhoI();
            obj.createComplianceBoundConstraintRhoD();
            obj.createVolumeConstraintWithBound();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,2,100);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.fHandle  = @(x) ones(size(x(1,:,:)));
            s.ndimf    = 1;
            s.mesh     = obj.mesh;
            aFun       = AnalyticalFunction(s);
            s.fun      = aFun.project('P1');
            s.mesh     = obj.mesh;
            s.type     = 'DensityAndBound';
            s.plotting = true;
            dens       = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilterRhoE(obj)
            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.75; % ERODED
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterRhoE = f;
        end

        function createFilterRhoI(obj)
            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.5; % INTERMEDIATE
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterRhoI = f;
        end

        function createFilterRhoD(obj)
            s.filterType = 'FilterAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.25; % DILATED
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterRhoD = f;
        end

        function createFilterAdjointRhoE(obj)
            s.filterType = 'FilterAdjointAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.75; % ERODED
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterGradientRhoE = f;
        end

        function createFilterAdointRhoI(obj)
            s.filterType = 'FilterAdjointAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.5; % INTERMEDIATE
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterGradientRhoI = f;
        end

        function createFilterAdjointRhoD(obj)
            s.filterType = 'FilterAdjointAndProject';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.filterStep = 'PDE';
            s.eta        = 0.25; % DILATED
            s.beta       = 1;
            f            = Filter.create(s);
            obj.filterGradientRhoD = f;
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            x = obj.designVariable.density.fun;       
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createLinearBoundFunction(obj)
            obj.linearBoundFunction = LinearBoundFunction();
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createComplianceBoundConstraintRhoE(obj)
            s.mesh                        = obj.mesh;
            s.filterDesignVariable        = obj.filterRhoE;
            s.filterGradient              = obj.filterGradientRhoE;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceWithBoundConstraint(s);
            obj.complianceRhoE = c;
        end

        function createComplianceBoundConstraintRhoI(obj)
            s.mesh                        = obj.mesh;
            s.filterDesignVariable        = obj.filterRhoI;
            s.filterGradient              = obj.filterGradientRhoI;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceWithBoundConstraint(s);
            obj.complianceRhoI = c;
        end

        function createComplianceBoundConstraintRhoD(obj)
            s.mesh                        = obj.mesh;
            s.filterDesignVariable        = obj.filterRhoD;
            s.filterGradient              = obj.filterGradientRhoD;
            s.complainceFromConstitutive  = obj.createComplianceFromConstiutive();
            s.material                    = obj.createMaterial();
            c = ComplianceWithBoundConstraint(s);
            obj.complianceRhoD = c;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraintWithBound(obj)
            s.mesh         = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v              = VolumeConstraintWithBound(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.linearBoundFunction;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*eye(n);
            M(n+1,n+1) = 1;
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.complianceRhoE;
            s.shapeFunctions{2} = obj.complianceRhoI;
            s.shapeFunctions{3} = obj.complianceRhoD;
            s.shapeFunctions{4} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 4;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
            s.ub             = [ones(obj.mesh.nnodes,1);1000];
            s.lb             = [zeros(obj.mesh.nnodes,1);-1000];
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

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