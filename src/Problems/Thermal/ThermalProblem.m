classdef ThermalProblem < handle
    
    properties (Access = public)
        uFun
        forces
    end

    properties (Access = private)
        boundaryConditions, bcApplier

        stiffness
        solverType, solverMode, solverCase

        problemSolver

        conductivity
        test
        trial
        source
    end

    properties (Access = protected)
        mesh 
        material  
    end

    methods (Access = public)

        function obj = ThermalProblem(cParams)
            obj.init(cParams);
            obj.createTemperatureFun(); 
            obj.createBCApplier();
            obj.createSolver();
        end

        function solve(obj, kappa)
            obj.computeStiffnessMatrix(kappa); % LHS
            obj.computeForces();          % RHS
            obj.computeTemperature();     % Solve PDE 
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.conductivity = cParams.conductivity;
            obj.source       = cParams.source;
            obj.solverType   = cParams.solverType;
            obj.solverMode   = cParams.solverMode;   % check what it means!
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;
            obj.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.trial = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function createTemperatureFun(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end

        function createSolver(obj)
            sS.type              = obj.solverCase;
            solver               = Solver.create(sS);
            s.solverType         = obj.solverType;
            s.solverMode         = obj.solverMode;
            s.solver             = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = obj.bcApplier;
            obj.problemSolver    = ProblemSolver(s);
        end

        function computeStiffnessMatrix(obj, kappa)
            s.test  = obj.test;
            s.trial = obj.trial;
            s.mesh  = obj.mesh;
            s.quadratureOrder = 2;
            s.function        = kappa; 
            s.type            = 'StiffnessMatrixWithFunction';
            lhs = LHSIntegrator.create(s);
            obj.stiffness = lhs.compute();
        end

        function computeForces(obj)
            s.type     = 'ShapeFunction';
            s.mesh     = obj.mesh;
            s.quadType = 2;
            RHSint = RHSIntegrator.create(s);
            rhs = RHSint.compute(obj.source, obj.test);
            obj.forces = rhs;
        end

        function u = computeTemperature(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.problemSolver.solve(s);           
            obj.uFun.setFValues(u);
        end

    end

end