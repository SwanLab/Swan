classdef GlobalPerimeterValidation < handle

    properties (Access = private)
        mesh
        designVariable
        perimeter
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = GlobalPerimeterValidation()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createPerimeter();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();

            saveas(gcf,'Paper/Validation/MonitoringGlobalPerimeterValidation.fig');
            obj.designVariable.fun.print('Paper/Validation/GlobalPerimeterValidationfValues');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(1,1,100,100);
        end

        function createDesignVariable(obj)
            s.type        = 'CircleInclusion';
            s.radius      = 0.15;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);

            %s.fun              = lsFun;
            s.fun              = LagrangianFunction.create(obj.mesh,1,'P1');
            s.fun.setFValues(1-heaviside(lsFun.fValues));
            s.mesh             = obj.mesh;
            s.type             = 'Density';%'LevelSet';
            s.plotting         = true;
            ls                 = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createPerimeter(obj)
            sF.mesh        = obj.mesh;
            sF.trial       = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.senseVector = ConstantFunction.create([0;1],obj.mesh);
            sF.ovAngleDeg  = 45;
            sF.tol         = 3e-4;
            filter         = FilterRegularizedDiamond(sF);

            h         = obj.mesh.computeMeanCellSize();
            s.mesh    = obj.mesh;
            s.epsilon = 8*h;
            s.value0 = 1;

            s.uMesh          = obj.createBaseDomain();
            s.filter         = filter;
            obj.perimeter = PerimeterFunctional(s);

        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.85;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.perimeter;
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

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            s.ub = 1;
            s.lb = 0;
            s.tauMax = 1000;
            s.tau = [];
            obj.primalUpdater = ProjectedGradient(s);%SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.etaNorm        = 0.1;
            s.etaNormMin     = 0.1;
            s.gJFlowRatio    = 0.5;
            s.etaMax         = 100;
            s.etaMaxMin      = 0.1;
            s.primalUpdater  = obj.primalUpdater;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end