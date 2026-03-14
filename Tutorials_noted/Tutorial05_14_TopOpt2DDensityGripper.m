classdef Tutorial05_14_TopOpt2DDensityGripper < handle

    properties (Access = private)
        filename
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        adjointProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        primalUpdater
        optimizer
        gJ
    end

    methods (Access = public)

        function obj = Tutorial05_14_TopOpt2DDensityGripper()
            obj.init();
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createAdjointProblem();
            obj.createNonSelfAdjCompliance();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'Gripping'; 
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            s.fun     = aFun.project('P1');
            s.mesh    = obj.mesh;
            s.type = 'Density';
            s.plotting = true;
            dens    = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

        function createMaterialInterpolator(obj)
            E0   = 1e-3;
            nu0  = 1/3;
            E1   = 1;
            nu1  = 1/3;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
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
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createAdjointProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditionsAdjoint();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.adjointProblem = fem;
        end

        function createNonSelfAdjCompliance(obj)
            s.mesh           = obj.mesh;
            s.filter         = obj.filter;
            s.material       = obj.createMaterial();
            s.stateProblem   = obj.physicalProblem;
            s.adjointProblem = obj.adjointProblem;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.4;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            s.tau    = [];
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 20;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 0.1;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filter.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh = obj.mesh;
            m = Material.create(s);
        end

        function bc = createBoundaryConditions(obj)
            isDir   = @(coor)  coor(:,2)>=0.49 & coor(:,2)<=0.51 & coor(:,1)<=0.1+1e-8 | coor(:,2)>=0.49 & coor(:,2)<=0.51 & coor(:,1)>=1-1e-8; % middle point of the gripping right part

            isPLTopRight      = @(coor)  (abs(coor(:,1)) >= 0.92 & coor(:,2) == 1 ); % top right part of the domain
            isPLBottomRight   = @(coor)  (abs(coor(:,1)) >= 0.92 & abs(coor(:,2)) <= 1e-8 );% bottom right part of the domain % not exactly 0 in the mesh
            isPLTopGripper    = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.6 ); % top part of the gripping zone
            isPLBottomGripper = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.4 ); % bottom part of the gripping zone

            sDir{1}.domain    = @(coor) isDir(coor); % fixed
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isPLTopRight(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -10; % downward force on the top right part

            sPL{2}.domain    = @(coor) isPLBottomRight(coor);
            sPL{2}.direction = 2;
            sPL{2}.value     = +10; % upward force on the bottom right part

            sPL{3}.domain    = @(coor) isPLTopGripper(coor);
            sPL{3}.direction = 2;
            sPL{3}.value     = +1; % upward force on the top part of the gripping zone (reaction)

            sPL{4}.domain    = @(coor) isPLBottomGripper(coor);
            sPL{4}.direction = 2;
            sPL{4}.value     = -1; % downward force on the bottom part of the gripping zone (reaction)

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

        function bc = createBoundaryConditionsAdjoint(obj)
            isDir   = @(coor)  coor(:,2)>=0.49 & coor(:,2)<=0.51 & coor(:,1)<=0.1+1e-8 | coor(:,2)>=0.49 & coor(:,2)<=0.51 & coor(:,1)>=1-1e-8;

            isPLTopRight      = @(coor)  (abs(coor(:,1)) >= 0.92 & coor(:,2) == 1 );
            isPLBottomRight   = @(coor)  (abs(coor(:,1)) >= 0.92 & abs(coor(:,2)) <= 1e-8 ); % not exactly 0 in the mesh
            isPLTopGripper    = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.6 );
            isPLBottomGripper = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.4 );

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isPLTopGripper(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = +2;

            sPL{2}.domain    = @(coor) isPLBottomGripper(coor);
            sPL{2}.direction = 2;
            sPL{2}.value     = -2;

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