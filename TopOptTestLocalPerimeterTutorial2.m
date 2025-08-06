classdef TopOptTestLocalPerimeterTutorial2 < handle

    properties (Access = private)
        mesh
        filterPerimeter
        designVariable
        volume
        perimeterConstraint
        cost
        constraint
        dualVariable
        optimizerGP
        locPer1
        locPer2
        locPer3
        locPer4
    end

    methods (Access = public)

        function obj = TopOptTestLocalPerimeterTutorial2()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterPerimeter();
            obj.createVolume();
            obj.createGlobalPerimeterConstraint();
            obj.createCostSimulations();
            obj.createConstraintSimulation1();
            obj.createDualVariableSimulation1();
            %obj.createOptimizerSimulation1();
            obj.createLocalPerimeterConstraints();
            obj.createConstraintSimulation2();
            obj.createDualVariableSimulation2();
            obj.createOptimizerSimulation2();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,1,75);
            x2      = linspace(0,1,75);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.mesh = Mesh.create(s);
        end

        function createDesignVariable(obj)
            sG.type    = 'Given';
            x1         = @(x) x(1,:,:);
            x2         = @(x) x(2,:,:);
            squareTL   = @(x) (x1(x)-0.25).^2+(x2(x)-0.75).^2-0.125^2;
            squareTR   = @(x) (x1(x)-0.75).^2+(x2(x)-0.75).^2-0.2^2;
            squareBL   = @(x) (x1(x)-0.25).^2+(x2(x)-0.25).^2-0.075^2;
            squareBR   = @(x) (x1(x)-0.75).^2+(x2(x)-0.25).^2-0.05^2;
            sG.fHandle = @(x) squareTL(x).*squareTR(x).*squareBL(x).*squareBR(x);
            g          = GeometricalFunction(sG);
            lsFun      = g.computeLevelSetFunction(obj.mesh);
            s.fun      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.fun.fValues = heaviside(lsFun.fValues); % Because it is an inclusion
            s.mesh     = obj.mesh;
            s.type     = 'Density';
            s.plotting = true;
            dens       = DesignVariable.create(s);
            obj.designVariable = dens;
        end

        function createFilterPerimeter(obj)
            s.filterType   = 'PDE';
            s.boundaryType = 'Neumann';
            s.metric       = 'Isotropy';
            s.mesh         = obj.mesh;
            s.trial        = LagrangianFunction.create(obj.mesh,1,'P1');
            f              = Filter.create(s);
            obj.filterPerimeter = f;
        end

        function createVolume(obj)
            s.mesh         = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            v              = VolumeFunctional(s);
            obj.volume     = v;
        end

        function createGlobalPerimeterConstraint(obj)
            s.mesh    = obj.mesh;
            s.filter  = obj.filterPerimeter;
            s.epsilon = 8*obj.mesh.computeMeanCellSize();
            s.minEpsilon = 8*obj.mesh.computeMeanCellSize();
            s.perimeterTargetAbs = 4*(0.125*pi); % Internal perimeter of 4 circular holes
            P         = PerimeterConstraint(s);
            obj.perimeterConstraint = P;
        end

        function createCostSimulations(obj)
            s.shapeFunctions{1} = obj.volume;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            LHS = LHSintegrator.create(s);
            M = LHS.compute;

            h = obj.mesh.computeMinCellSize();
            M = h^2*eye(size(M));
        end

        function createConstraintSimulation1(obj)
            s.shapeFunctions{1} = obj.perimeterConstraint;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariableSimulation1(obj)
            s.nConstraints   = 1;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizerSimulation1(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 600;
            s.tolerance      = 1e-8;
            s.constraintCase = {'EQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 1;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizerGP = opt;
        end

        function createLocalPerimeterConstraints(obj)
            x1 = @(x) x(1,:,:);
            x2 = @(x) x(2,:,:);
            s.mesh    = obj.mesh;
            s.filter  = obj.filterPerimeter;
            s.epsilon = 8*obj.mesh.computeMeanCellSize();
            s.minEpsilon = 4*obj.mesh.computeMeanCellSize();
            s.perimeterTargetAbs = 0.25*pi;
            s.subDomainHandle = @(x) max(abs(x1(x)-0.25),abs(x2(x)-0.75))/0.5 - 0.5;
            P         = LocalPerimeterConstraint(s);
            obj.locPer1 = P;

            s.perimeterTargetAbs = 0.25*pi;
            s.subDomainHandle = @(x) max(abs(x1(x)-0.75),abs(x2(x)-0.75))/0.5 - 0.5;
            P         = LocalPerimeterConstraint(s);
            obj.locPer2 = P;

            s.perimeterTargetAbs = 0.25*pi;
            s.subDomainHandle = @(x) max(abs(x1(x)-0.25),abs(x2(x)-0.25))/0.5 - 0.5;
            P         = LocalPerimeterConstraint(s);
            obj.locPer3 = P;

            s.perimeterTargetAbs = 0.25*pi;
            s.subDomainHandle = @(x) max(abs(x1(x)-0.75),abs(x2(x)-0.25))/0.5 - 0.5;
            P         = LocalPerimeterConstraint(s);
            obj.locPer4 = P;
        end

        function createConstraintSimulation2(obj)
            s.shapeFunctions{1} = obj.locPer1;
            s.shapeFunctions{2} = obj.locPer2;
            s.shapeFunctions{3} = obj.locPer3;
            s.shapeFunctions{4} = obj.locPer4;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariableSimulation2(obj)
            s.nConstraints   = 4;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizerSimulation2(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 600;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 1;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizerGP = opt;
        end

    end
end