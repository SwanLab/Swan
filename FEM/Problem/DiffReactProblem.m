classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
    end
    
    properties (Access = private)
        dim
        mesh
        solver
        epsilon
        LHStype
        problemLHS
        problemData
        boundaryConditions
    end

    methods (Access = public)
        
        function obj = DiffReactProblem(cParams)
            obj.init(cParams);
            obj.computeDimensions();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.createProblemLHS();
        end
        
        function solve(obj)
            obj.epsilon = .1857;
            bc  = obj.boundaryConditions;
            M = obj.problemLHS.K;
            rhs = M*ones(size(M,1), 1);
            RHS = bc.fullToReducedVector(rhs);
            LHS = obj.computeLHS();
            x = obj.solver.solve(LHS,RHS);
            obj.variables.x = bc.reducedToFullVector(x);
        end

        function computeVariables(obj,rhs)
            bc  = obj.boundaryConditions;
            RHS = bc.fullToReducedVector(rhs);
            LHS = obj.computeLHS(obj.epsilon);
            x = obj.solver.solve(LHS,RHS);
            obj.variables.x = bc.reducedToFullVector(x);
        end
        
        function LHS = computeLHS(obj, epsilon)
            obj.epsilon = epsilon;
            lhs = obj.problemLHS.compute(epsilon);
            LHS = obj.boundaryConditions.fullToReducedMatrix(lhs);
        end
       
        function print(obj,filename)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            s.quad = quad;
            s.mesh = obj.mesh;
            s.iter = 0;
            s.fields    = obj.variables.x;
            s.ptype     = 'DIFF-REACT';
            s.ndim      = 2;
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

        function computeDimensions(obj)
            s.name = 'x';
            s.mesh = obj.mesh;
            dims   = DimensionScalar(s);
            obj.dim = dims;
        end

        function createBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.scale        = obj.problemData.scale;
            s.bc.dirichlet = [];
            s.bc.pointload = [];
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end
        
        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function createProblemLHS(obj)
            s.type         = obj.LHStype;
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = [];
            obj.problemLHS = LHSintegrator.create(s);
        end
    
    end

end