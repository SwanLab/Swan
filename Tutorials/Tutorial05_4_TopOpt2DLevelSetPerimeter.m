classdef Tutorial05_4_TopOpt2DLevelSetPerimeter < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        perimeter
        volume
        cost
        constraint
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = Tutorial05_4_TopOpt2DLevelSetPerimeter()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createPerimeter();
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
            obj.mesh = TriangleMesh(1,1,50,50);
        end

        function createDesignVariable(obj)
            s.type        = 'SquareInclusion';
            s.length      = 0.5;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            s.fun         = lsFun;
            s.mesh        = obj.mesh;
            s.type        = 'LevelSet';
            s.plotting    = true;
            ls            = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            u = 65;
            alpha = 90;
            CAnisotropic = [tand(u),0;0,1/tand(u)];
            R = [cosd(alpha),-sind(alpha)
                sind(alpha), cosd(alpha)];
            CGlobal = R*CAnisotropic*R';

            s.filterType   = 'PDE';
            s.boundaryType = 'Neumann';
            s.metric       = 'Anisotropy';
            s.A            = ConstantFunction.create(CGlobal,obj.mesh);
            s.mesh         = obj.mesh;
            s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            f              = Filter.create(s);
            obj.filter     = f;
        end

        function createPerimeter(obj)
            eOverhmin     = 10;
            epsilon       = eOverhmin*obj.mesh.computeMeanCellSize();
            s.mesh        = obj.mesh;
            s.filter      = obj.filter;
            s.epsilon     = epsilon;
            s.value0      = 4; % external P
            s.uMesh       = obj.createBaseDomain();
            P             = PerimeterFunctional(s);
            obj.perimeter = P;
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
            s.GIFname = 'per';
            s.designVariable = obj.designVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.volumeTarget   = 0.85;
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 5;
            s.etaMax         = 1;
            s.etaMaxMin      = 0.01;
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