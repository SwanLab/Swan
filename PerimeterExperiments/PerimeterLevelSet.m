classdef PerimeterLevelSet < handle

    properties (Access = private)
        mesh
        designVariable
        volume
        filter
        perimeter
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = PerimeterLevelSet()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createVolumeConstraint();
            obj.createFilter();
            obj.createPerimeter();
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
            obj.mesh = TriangleMesh(1,1,200,200);
        end

        function createDesignVariable(obj)
            s.type = 'SquareInclusion';
            s.length = 0.35;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g = GeometricalFunction(s);
            ls = g.computeLevelSetFunction(obj.mesh);
            
            sD.fun      = LagrangianFunction.create(obj.mesh,1,'P1');
            sD.fun.setFValues(ls.fValues);
            sD.mesh     = obj.mesh;
            sD.type     = 'LevelSet';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
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
            s.volumeTarget = 0.85;
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createFilter(obj)
            s.mesh  = obj.mesh;
            s.alpha = 4;
            s.beta  = 0;
            s.theta = 90;
            s.trial  = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.filter  = FilterPDE(s);
        end

        function createPerimeter(obj)
            h         = obj.mesh.computeMinCellSize();
            s.mesh    = obj.mesh;
            s.uMesh   = obj.createBaseDomain();
            s.filter  = obj.filter;
            s.epsilon = 10*h;
            s.value0  = 1;
            obj.perimeter = PerimeterFunctional(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.perimeter;
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
            s.shapeFunctions{1} = obj.volume;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 3000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.etaNorm        = 0.01;
            s.etaNormMin = 0.01;
            s.etaMax = 10;
            s.etaMaxMin = 0.1;
            s.gif = false;
            s.gifName = [];
            s.printing = true;
            s.printName = 'Results/LevelSetSegment';
            s.gJFlowRatio    = 2;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end