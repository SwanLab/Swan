classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
        x
    end
    
    properties (Access = private)
        mesh
        solver
        epsilon
        LHStype
        problemLHS
        problemData
        % boundaryConditions
    end

    methods (Access = public)
        
        function obj = DiffReactProblem(cParams)
            obj.init(cParams);
            % obj.createBoundaryConditions();
            obj.createSolver();
            obj.createProblemLHS();
        end

        function computeVariables(obj,rhs)
            % bc  = obj.boundaryConditions;
            % RHS = bc.fullToReducedVector(rhs);
            RHS = rhs;
            LHS = obj.computeLHS(obj.epsilon);
            x = obj.solver.solve(LHS,RHS);
            % obj.variables.x = bc.reducedToFullVector(x);
            obj.variables.x = x;
            a.mesh = obj.mesh;
            a.fValues = obj.variables.x;
            a.order = 'P1';
            obj.x = LagrangianFunction(a);
        end
        
        function LHS = computeLHS(obj, epsilon)
            obj.epsilon = epsilon;
            lhs = obj.problemLHS.compute(epsilon);
            % LHS = obj.boundaryConditions.fullToReducedMatrix(lhs);
            LHS = lhs;
        end
       
        function print(obj,filename)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            s.quad = quad;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.variables.x;
            s.ptype     = 'DIFF-REACT';
            s.ndim      = 3;
            s.pdim      = obj.problemData.pdim;
            s.type      = 'ScalarNodal';
            fPrinter = FemPrinter(s);
            fPrinter.print(filename);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh              = cParams.mesh;
            obj.LHStype           = cParams.LHStype;
            obj.problemData.pdim  = '1D';
            obj.problemData.scale = cParams.scale;
        end

        function createBoundaryConditions(obj)
            s.mesh  = obj.mesh;
            s.scale = obj.problemData.scale;
            s.ndofs = obj.mesh.nnodes;
            s.bc{1}.dirichlet = [];
            s.bc{1}.pointload = [];
            s.bc{1}.ndimf     = [];
            s.bc{1}.ndofs     = [];
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end
        
        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end

        function createProblemLHS(obj)
            s.type  = obj.LHStype;
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.quadratureOrder = 2;
            obj.problemLHS = LHSintegrator.create(s);
        end
    
    end

end