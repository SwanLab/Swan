classdef ProblemConstructorReducedDisp < handle

    properties (Access = private)
        boundaryConditions
        solver
    end
    
    methods (Access = public)
        
        function obj = ProblemConstructorReducedDisp(cParams)
            obj.init(cParams)
        end
        
        function [u,L] = solve(obj,cParams)
            K = cParams.stiffness;
            F = cParams.forces; %%% This is not inside BC???
            [LHS, RHS] = obj.computeMatrices(K,F);
            sol        = obj.solver.solve(LHS, RHS);
            [u, L]     = obj.cleanSolution(sol,K);
        end
   end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solver             = cParams.solver;
        end
        
        function [LHS,RHS] = computeMatrices(obj,K,F)
            LHS = obj.assembleLHS(K);
            RHS = obj.assembleRHS(K,F);
        end
        
        function [LHS] = assembleLHS(K)
            bcs = obj.boundaryConditions;
            dofs = 1:size(K);
            free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
            LHS = K(free_dofs, free_dofs);
        end
        
        function [RHS] = assembleRHS(K,F)
            bcs = obj.boundaryConditions;
            dofs = 1:size(K);
            free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
            RHS = F(free_dofs);
        end
        
        function [u,L] = cleanSolution(obj,sol,K)
            dofs = 1:size(K);
            free_dofs = setdiff(dofs, bcs.dirichlet_dofs);
            u = zeros(size(K,1), 1);
            u(free_dofs) = sol;
            u(bcs.dirichlet_dofs) = bcs.dirichlet_vals;
            L = [];
        end
    end
    
end