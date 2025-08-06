classdef TopOptTestLocalPerimeterTutorial < handle

    properties (Access = private)
        mesh
        filterPerimeter
        designVariable
        volume
        perimeterConstraintL
        perimeterConstraintR
        cost
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = TopOptTestLocalPerimeterTutorial()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilterPerimeter();
            obj.createVolume();
            obj.createPerimeterConstraintLeft();
            obj.createPerimeterConstraintRight();
            obj.createCostSimulations();
            obj.createConstraintSimulation();
            obj.createDualVariableSimulation();
            obj.createOptimizerSimulation();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            %UnitMesh better
            x1      = linspace(0,1,151);
            x2      = linspace(0,1,151);
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
            squareTL   = @(x) max(abs(x1(x)-0.25),abs(x2(x)-0.75))/0.25 - 0.5;
            squareTR   = @(x) max(abs(x1(x)-0.75),abs(x2(x)-0.75))/0.25 - 0.5;
            squareBL   = @(x) max(abs(x1(x)-0.25),abs(x2(x)-0.25))/0.25 - 0.5;
            squareBR   = @(x) max(abs(x1(x)-0.75),abs(x2(x)-0.25))/0.25 - 0.5;
            sG.fHandle = @(x) squareTL(x).*squareTR(x).*squareBL(x).*squareBR(x);
            g          = GeometricalFunction(sG);
            lsFun      = g.computeLevelSetFunction(obj.mesh);
            s.fun      = LagrangianFunction.create(obj.mesh,1,'P1');
            s.fun.setFValues(heaviside(lsFun.fValues)); % Because it is an inclusion
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

        function createPerimeterConstraintLeft(obj)
            s.mesh       = obj.mesh;
            s.uMesh      = obj.computeSubDomain(0.25,0.5,0.5,1);
            s.filter     = obj.filterPerimeter;
            s.epsilon    = 6*obj.mesh.computeMeanCellSize();
            s.minEpsilon = 2*obj.mesh.computeMeanCellSize();
            s.value0     = 1;
            s.target     = 2*(0.125*pi); % Internal perimeter of 2 circular holes
            P            = PerimeterConstraint(s);
            obj.perimeterConstraintL = P;
        end

        function createPerimeterConstraintRight(obj)
            s.mesh       = obj.mesh;
            s.uMesh      = obj.computeSubDomain(0.75,0.5,0.5,1);
            s.filter     = obj.filterPerimeter;
            s.epsilon    = 6*obj.mesh.computeMeanCellSize();
            s.minEpsilon = 2*obj.mesh.computeMeanCellSize();
            s.value0     = 1;
            s.target     = 2*(0.25*pi); % Internal perimeter of 2 circular holes
            P            = PerimeterConstraint(s);
            obj.perimeterConstraintR = P;
        end

        function uMesh = computeSubDomain(obj,x0,y0,dx,dy)
            sG.type            = 'Rectangle';
            sG.xCoorCenter     = x0;
            sG.yCoorCenter     = y0;
            sG.xSide           = dx;
            sG.ySide           = dy;
            g                  = GeometricalFunction(sG);
            lsFun              = g.computeLevelSetFunction(obj.mesh);
            sUm.backgroundMesh = obj.mesh;
            sUm.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsFun.fValues);
        end

        function createCostSimulations(obj)
            s.shapeFunctions{1} = obj.volume;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraintSimulation(obj)
            s.shapeFunctions{1} = obj.perimeterConstraintL;
            s.shapeFunctions{2} = obj.perimeterConstraintR;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createDualVariableSimulation(obj)
            s.nConstraints   = 2;
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function p = createPrimalUpdater(obj)
            s.ub     = 1;
            s.lb     = 0;
            s.tauMax = 1000;
            p        = ProjectedGradient(s);
        end

        function createOptimizerSimulation(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 600;
            s.tolerance      = 1e-8;
            s.constraintCase = {'INEQUALITY','INEQUALITY'};
            s.primal         = 'PROJECTED GRADIENT';
            s.ub             = 1;
            s.lb             = 0;
            s.etaNorm        = 0.01;
            s.gJFlowRatio    = 1;
            s.primalUpdater  = obj.createPrimalUpdater();
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end
    end
end
