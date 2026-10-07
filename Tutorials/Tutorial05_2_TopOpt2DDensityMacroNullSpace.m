classdef Tutorial05_2_TopOpt2DDensityMacroNullSpace < handle

    properties (Access = private)
        mesh
        filter
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    properties (Access = public)
        designVariable
        rhoSmooth
    end

    methods (Access = public)

        function hist = getCostHistory(obj)
            hist = obj.optimizer.costHistory;
        end

        function hist = getComplianceHistory(obj)
            hist = obj.optimizer.complianceHistory;
        end

        function hist = getVolumeConstraintHistory(obj)
            hist = obj.optimizer.volumeConstraintHistory;
        end

        function rho = getRhoValues(obj)
            rho = obj.designVariable.fun.fValues;
        end

        function rhoSmooth = getRhoSmoothValues(obj)
            rhoSmooth = obj.rhoSmooth.fValues;
        end

        function meshData = getMeshData(obj)

            meshData.coord  = obj.mesh.coord;
            meshData.connec = obj.mesh.connec;
            meshData.ndim   = obj.mesh.ndim;
            meshData.nnodes = obj.mesh.nnodes;

        end

        function obj = Tutorial05_2_TopOpt2DDensityMacroNullSpace()

            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createComplianceFromConstiutive();
            obj.createCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            obj.rhoSmooth = ...
                obj.filter.compute( ...
                obj.designVariable.fun,2);

        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)

            obj.mesh = TriangleMesh(2,1,150,75);

        end

        function createDesignVariable(obj)

            s.fHandle = @(x) ...
                0.998*ones(size(x(1,:,:)));

            s.ndimf = 1;
            s.mesh  = obj.mesh;

            aFun = AnalyticalFunction(s);

            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;

            obj.designVariable = ...
                DesignVariable.create(sD);

        end

        function createFilter(obj)

            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = ...
                LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            obj.filter = Filter.create(s);

        end

        function createMaterialInterpolator(obj)

            E0 = 1e-3;
            nu0 = 0.3;

            ndim = obj.mesh.ndim;

            matA.shear = ...
                IsotropicElasticMaterial.computeMuFromYoungAndPoisson( ...
                E0,nu0);

            matA.bulk = ...
                IsotropicElasticMaterial.computeKappaFromYoungAndPoisson( ...
                E0,nu0,ndim);

            E1 = 1;
            nu1 = 0.3;

            matB.shear = ...
                IsotropicElasticMaterial.computeMuFromYoungAndPoisson( ...
                E1,nu1);

            matB.bulk = ...
                IsotropicElasticMaterial.computeKappaFromYoungAndPoisson( ...
                E1,nu1,ndim);

            s.interpolation = 'SIMPALL';
            s.dim           = '2D';
            s.matA          = matA;
            s.matB          = matB;

            obj.materialInterpolator = ...
                MaterialInterpolator.create(s);

        end

        function m = createMaterial(obj)

            f = obj.designVariable.fun;

            s.type = 'DensityBased';

            s.density = f;

            s.materialInterpolator = ...
                obj.materialInterpolator;

            s.dim  = '2D';
            s.mesh = obj.mesh;

            m = Material.create(s);

        end

        function createElasticProblem(obj)

            s.mesh     = obj.mesh;
            s.scale    = 'MACRO';
            s.material = obj.createMaterial();
            s.dim      = '2D';

            s.boundaryConditions = ...
                obj.createBoundaryConditions();

            s.interpolationType = 'LINEAR';
            s.solverType        = 'REDUCED';
            s.solverMode        = 'DISP';
            s.solverCase        = DirectSolver();

            obj.physicalProblem = ...
                ElasticProblem(s);

        end

        function c = createComplianceFromConstiutive(obj)

            s.mesh = obj.mesh;

            s.stateProblem = ...
                obj.physicalProblem;

            c = ...
                ComplianceFromConstitutiveTensor(s);

        end

        function createCompliance(obj)

            s.mesh   = obj.mesh;
            s.filter = obj.filter;

            s.complainceFromConstitutive = ...
                obj.createComplianceFromConstiutive();

            s.material = ...
                obj.createMaterial();

            obj.compliance = ...
                ComplianceFunctional(s);

        end

        function uMesh = createBaseDomain(obj)

            levelSet = ...
                -ones(obj.mesh.nnodes,1);

            s.backgroundMesh = obj.mesh;

            s.boundaryMesh = ...
                obj.mesh.createBoundaryMesh();

            uMesh = UnfittedMesh(s);

            uMesh.compute(levelSet);

        end

        function createVolumeConstraint(obj)

            s.mesh   = obj.mesh;
            s.filter = obj.filter;

            s.test = ...
                LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            s.volumeTarget = 0.4;

            s.uMesh = ...
                obj.createBaseDomain();

            obj.volume = ...
                VolumeConstraint(s);

        end

        function createCost(obj)

            s.shapeFunctions{1} = ...
                obj.compliance;

            s.weights = 1;

            s.Msmooth = ...
                obj.createMassMatrix();

            obj.cost = Cost(s);

        end

        function M = createMassMatrix(obj)

            test = ...
                LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            trial = ...
                LagrangianFunction.create( ...
                obj.mesh,1,'P1');

            M = ...
                IntegrateLHS( ...
                @(u,v) DP(v,u), ...
                test, ...
                trial, ...
                obj.mesh, ...
                'Domain');

            M = diag(sum(M,1));

        end

        function createConstraint(obj)

            s.shapeFunctions{1} = ...
                obj.volume;

            s.Msmooth = ...
                obj.createMassMatrix();

            obj.constraint = ...
                Constraint(s);

        end

        function createPrimalUpdater(obj)

            s.ub = 0.998;
            s.lb = 1e-6;

            s.tauMax = 100;
            s.tau    = [];

            obj.primalUpdater = ...
                ProjectedGradient(s);

        end

        function createOptimizer(obj)

            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;

            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};

            s.primal        = ...
                'PROJECTED GRADIENT';

            s.etaNorm       = 0.01;
            s.gJFlowRatio   = 2;
            s.primalUpdater = ...
                obj.primalUpdater;

            s.gif       = false;
            s.gifName   = [];
            s.printing  = false;
            s.printName = [];

            opt = OptimizerNullSpace(s);

            opt.solveProblem();

            obj.optimizer = opt;

        end

        function bc = createBoundaryConditions(obj)

            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDir = @(coor) ...
                abs(coor(:,1)-xMin) < 1e-12;

            isForce = @(coor) ...
                abs(coor(:,1)-xMax) < 1e-12 & ...
                coor(:,2) >= 0.4*yMax & ...
                coor(:,2) <= 0.6*yMax;

            sDir{1}.domain = ...
                @(coor) isDir(coor);

            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain = ...
                @(coor) isForce(coor);

            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dirichletFun = ...
                    [dirichletFun, ...
                    DirichletCondition( ...
                    obj.mesh,sDir{i})];

            end

            pointloadFun = [];

            for i = 1:numel(sPL)

                pointloadFun = ...
                    [pointloadFun, ...
                    TractionLoad( ...
                    obj.mesh, ...
                    sPL{i}, ...
                    'DIRAC')];

            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = pointloadFun;
            s.periodicFun  = [];
            s.mesh         = obj.mesh;

            bc = BoundaryConditions(s);

        end

    end

end