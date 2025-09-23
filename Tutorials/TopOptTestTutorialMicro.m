classdef TopOptTestTutorialMicro < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        materialInterpolator
        physicalProblem
        ChomogAlphaBeta
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestTutorialMicro()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createMicroAlphaBeta();
            obj.createVolumeConstraint();
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
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            sG.type            = 'CircleInclusion';
            sG.xCoorCenter     = 0.5;
            sG.yCoorCenter     = 0.5;
            sG.radius          = 0.25;
            g                  = GeometricalFunction(sG);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sF.fValues         = 1-heaviside(lsFun.fValues);
            sF.mesh            = obj.mesh;
            sF.order           = 'P1';
            s.fun              = LagrangianFunction(sF);
            s.mesh             = obj.mesh;
            s.type             = 'Density';
            s.plotting         = true;
            dens               = DesignVariable.create(s);
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
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = f{1}.project('P1');            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            s.solverCase = 'DIRECT';
            fem = ElasticProblemMicro(s);
            obj.physicalProblem = fem;
        end

        function createMicroAlphaBeta(obj)
            s.mesh              = obj.mesh;
            s.filter            = obj.filter;
            s.material          = obj.createMaterial();
            s.stateProblem      = obj.physicalProblem;
            s.alpha             = [0, 1; 1 ,0];
            s.beta              = [0, 1; 1 ,0];
            obj.ChomogAlphaBeta = MicroAlphaBetaFunctional(s);
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
            s.volumeTarget = 0.5;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.ChomogAlphaBeta;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            test  = LagrangianFunction.create(obj.mesh,1,'P1');
            trial = LagrangianFunction.create(obj.mesh,1,'P1'); 
            M = IntegrateLHS(@(u,v) DP(v,u),test,trial,obj.mesh);  
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
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
            s.maxIter        = 300;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;
            s.gJFlowRatio    = 0.2;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
            xMax = max(obj.mesh.coord(:,1));
            yMax = max(obj.mesh.coord(:,2));

            isDirBL = @(coor)  abs(coor(:,1))==0 & abs(coor(:,2))==0;
            isDirTL = @(coor)  abs(coor(:,1))==0 & abs(coor(:,2))==yMax;
            isDirTR = @(coor)  abs(coor(:,1))==xMax & abs(coor(:,2))==yMax;
            isDirBR = @(coor)  abs(coor(:,1))==xMax & abs(coor(:,2))==0;
            isDir   = @(coor) isDirBL(coor) | isDirTL(coor) | isDirTR(coor) | isDirBR(coor);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            isLeft  = @(coor) abs(coor(:,1))==0 & ~(isDirBL(coor)|isDirTL(coor));
            isRight = @(coor) abs(coor(:,1))==xMax & ~(isDirBR(coor)|isDirTR(coor));

            sPer{1}.leader   = @(coor) isLeft(coor);
            sPer{1}.follower = @(coor) isRight(coor);

            isBottom = @(coor) abs(coor(:,2))==0 & ~(isDirBL(coor)|isDirBR(coor));
            isTop    = @(coor) abs(coor(:,2))==yMax & ~(isDirTL(coor)|isDirTR(coor));

            sPer{2}.leader   = @(coor) isBottom(coor);
            sPer{2}.follower = @(coor) isTop(coor);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];

            periodicFun = [];
            for i = 1:numel(sPer)
                per = PeriodicCondition(obj.mesh,sPer{i});
                periodicFun = [periodicFun, per];
            end
            s.periodicFun  = PeriodicCondition;
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end