classdef DiffReactProblem < handle
    
    properties (GetAccess = public, SetAccess = protected)
        variables
    end
    
    properties (Access = private)
        mesh
        field
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
            obj.createField();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.createProblemLHS();
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
            s.ndim      = obj.mesh.ndim;
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

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'LINEAR';
            obj.field = Field(s);
        end

        function createBoundaryConditions(obj)
            s.dim   = obj.field.dim;
            s.mesh  = obj.mesh;
            s.scale = obj.problemData.scale;
            s.ndofs = obj.field.dim.ndofs;
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
            s.type = obj.LHStype;
            s.mesh = obj.mesh;
            obj.problemLHS = LHSintegrator.create(s);
        end
    
    end

end