classdef ConstraintSolverMonolitic < ConstraintSolverFactory

    properties (Access = private)
        solver
        bc
        K
        RHS
        sizeK
    end

    methods (Access = public)
        function obj = ConstraintSolverMonolitic(cParams)
            obj.init(cParams);
        end

        function [u, L] = solveSystem(obj, LHSMatrix, RHSMatrix, nConstraints)
            lhs = obj.createGeneralMatrix(LHSMatrix, nConstraints);
            sol = obj.solver.solve(lhs, RHSMatrix);
            u   = sol(1:obj.sizeK, 1);
            L   = -sol(obj.sizeK+1:end, 1);    
        end

    end

    methods (Access = private)
        function init(obj, cParams)
            obj.solver = cParams.solver;
            obj.bc     = cParams.bc; 
            obj.K      = cParams.LHS;
            obj.sizeK  = size(obj.K, 1);
            obj.RHS    = cParams.RHS;
        end

        function fullLHS = createGeneralMatrix(obj, Ct, nConstraints)
            C       = Ct';
            nC      = nConstraints;
            Z       = zeros(nC);
            Km      = obj.K;
            fullLHS = [Km C; C' Z];
        end

    end

end