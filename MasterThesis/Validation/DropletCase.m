classdef DropletCase < handle

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

        function obj = DropletCase()
            tic;
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
            disp(toc);
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            x1      = linspace(0,1,50);
            x2      = linspace(0,1,50);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            s.type        = 'CircleInclusion';
            s.radius      = 0.2;
            s.xCoorCenter = 0.5;
            s.yCoorCenter = 0.5;
            g             = GeometricalFunction(s);
            lsFun         = g.computeLevelSetFunction(obj.mesh);
            sD.fValues    = 1-heaviside(lsFun.fValues);
            sD.mesh       = obj.mesh;
            sD.order      = 'P1';
            s.fun         = LagrangianFunction(sD);
            s.mesh        = obj.mesh;
            s.type        = 'Density';
            s.plotting    = true;
            ls            = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.trial    = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh     = obj.mesh;
            s.type     = 'Droplet';
            s.theta    = 90;
            s.alpha    = 1.155;
            s.epsilon  = 1*obj.mesh.computeMeanCellSize();
            f          = NonLinearFilter.create(s);
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
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.6;
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
%             s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.mesh  = obj.mesh;
%             s.type  = 'MassMatrix';
%             LHS = LHSintegrator.create(s);
%             M = LHS.compute;
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

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.ub             = 1;
            s.lb             = 0;
            s.primal         = 'PROJECTED GRADIENT';
            s.etaNorm        = 0.005;
            s.gJFlowRatio    = 15;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end