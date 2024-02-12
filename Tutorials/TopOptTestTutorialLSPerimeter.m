classdef TopOptTestTutorialLSPerimeter < handle

    properties (Access = private)
        mesh
        filter
        designVariable
        perimeter
        volume
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestTutorialLSPerimeter()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();            
            obj.createFilter();
            obj.createPerimeter();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createOptimizer();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh(s);
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
            ls            = DesignVariable.create(s);   
            obj.designVariable = ls;
        end

        function createFilter(obj)
            u              = 65;
            alpha          = 90;
            s.filterType   = 'PDE';
            s.boundaryType = 'Neumann';
            s.metric       = 'Anisotropy';
            s.CAnisotropic = [tand(u),0;0,1/tand(u)];
            s.aniAlphaDeg  = alpha;
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end   

        function createPerimeter(obj)
            eOverhmin     = 10;
            epsilon       = eOverhmin*obj.mesh.computeMeanCellSize();
            s.mesh        = obj.mesh;
            s.filter      = obj.filter;
            s.epsilon     = epsilon;
            s.value0      = 4; % external P
            P             = PerimeterFunctional(s);
            obj.perimeter = P;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filter;
            s.volumeTarget = 0.85;
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.perimeter;
            s.weights           = 1;
            obj.cost            = Cost(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
            obj.constraint      = Constraint(s);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 40;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.volumeTarget   = 0.85;
            s.primal         = 'SLERP';
            s.aJmax          = 12;
            s.aGmax          = 3.5e5;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end