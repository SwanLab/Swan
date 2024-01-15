classdef ProblemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        type
        stiffness
        forces
        boundaryConditions
        BCApplier
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ProblemSolver(cParams)
            obj.init(cParams)
        end

        function u = solve(obj)
            [LHS, RHS] = obj.computeMatrices();
            sol        = obj.solveSystem(LHS, RHS);
            [u, L]     = obj.cleanupSolution(sol);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type               = cParams.type;
            obj.stiffness          = cParams.stiffness;
            obj.forces             = cParams.forces;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.BCApplier          = cParams.BCApplier;
        end

        function [LHS, RHS] = computeMatrices(obj)
            LHS  = obj.assembleLHS();
            RHS  = obj.assembleRHS();
        end

        function sol = solveSystem(obj, LHS, RHS)
            sol = LHS\RHS;
        end

        function [u, L] = cleanupSolution(obj,sol)
            bcapp = obj.BCApplier;
            switch obj.type
                case 'MONOLITHIC'
                    nDisp = size(obj.stiffness,1);
                    u = sol(1:nDisp, :);
                    L = sol( (nDisp+1):end, : );
                    
                case 'REDUCED'
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    u = zeros(size(obj.stiffness,1), 1);
                    u(free_dofs) = sol;
                    u(bcapp.dirichlet_dofs) = bcapp.dirichlet_vals;
                    L = [];
            end

            % maybe return a p1function or whatever
        end

        function LHS = assembleLHS(obj)
            bcapp = obj.BCApplier;
            switch obj.type
                case 'MONOLITHIC'
                    Ct = bcapp.computeLinearConditionsMatrix();
                    C   = Ct';
                    nC  = size(Ct,1);
                    Z   = zeros(nC);
                    Km  = obj.stiffness;
                    LHS = [Km C; C' Z];
                    
                case 'REDUCED'
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    LHS = obj.stiffness(free_dofs, free_dofs);

            end
        end

        function RHS = assembleRHS(obj)
            bcapp = obj.BCApplier;
            switch obj.type
                case 'MONOLITHIC'
                    dir_vals = bcapp.dirichlet_vals;
                    nCases = size(obj.forces,2);
                    Ct = repmat(dir_vals, [1 nCases]);
                    RHS = [obj.forces; Ct];
                    
                case 'REDUCED'
                    dofs = 1:size(obj.stiffness);
                    free_dofs = setdiff(dofs, bcapp.dirichlet_dofs);
                    RHS = obj.forces(free_dofs);

            end
        end

        
    end
    
end