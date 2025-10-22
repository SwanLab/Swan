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
            obj.computeForces();  % Case where RHS is fixed. If updated, call update function before 'solve'
        end

        function updateConductivity(obj,kappa)
            obj.computeStiffnessMatrix(kappa); % LHS
        end

        function updateRHSWithMass(obj, mass) 
            obj.forces = IntegrateRHS(@(v) mass.* DP(obj.source,v), obj.test, obj.mesh, 2);
        end

        function solve(obj)
            obj.computeTemperature();          % Solve PDE 
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
            obj.stiffness =IntegrateLHS(@(u,v) kappa.*DP(Grad(u),Grad(v)),obj.test,obj.trial,obj.mesh,'Domain',2);
        end

        function computeForces(obj)
            obj.forces = IntegrateRHS(@(v) DP(obj.source,v), obj.test, obj.mesh, 2);
        end

        function u = computeTemperature(obj)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces;
            [u,~]       = obj.problemSolver.solve(s);           
            obj.uFun.setFValues(u);
        end

    end

end
