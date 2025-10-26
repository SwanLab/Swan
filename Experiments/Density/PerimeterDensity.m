classdef PerimeterDensity < handle

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

        function obj = PerimeterDensity()
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
            sD.fun.setFValues(1-heaviside(ls.fValues));
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
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
            u = 70;
            alpha = 90;
            CAnisotropic = [tand(u),0;0,1/tand(u)];
            R = [cosd(alpha),-sind(alpha)
                sind(alpha), cosd(alpha)];
            CGlobal = R*CAnisotropic*R';

            sF.A     = ConstantFunction.create(CGlobal,obj.mesh);
            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh = obj.mesh;
            sF.filterType = 'PDE';
            sF.boundaryType = 'Neumann';
            sF.metric = 'Anisotropy';
            obj.filter = Filter.create(sF);
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
            s.maxIter        = 1200;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.etaNorm        = 0.1;
            s.gif = false;
            s.gifName = [];
            s.printing = true;
            s.printName = 'Results/DensityEllipseVertical';
            s.gJFlowRatio    = 0.2;
            s.primalUpdater  = obj.primalUpdater;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end