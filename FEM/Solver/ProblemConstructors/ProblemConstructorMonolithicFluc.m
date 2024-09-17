classdef ProblemConstructorMonolithicFluc < ProblemSolver

    properties (Access = private)
        bcApplier
    end

    methods (Access = public)
        
        function obj = ProblemConstructorMonolithicFluc(cParams)
            obj@ProblemSolver(cParams);
            obj.bcApplier          = cParams.BCApplier;
        end

        function [u,L] = solve(obj,s)
            [LHS, RHS] = obj.computeMatrices(s);
            sol        = obj.solver.solve(LHS, RHS);
            [u, L]     = obj.cleanSolution(sol,s);
        end
   end
    
    methods (Access = private)
        
        function [LHS,RHS] = computeMatrices(obj,s)
            LHS = obj.assembleLHS(s);
            RHS = obj.assembleRHS(s);
        end
        
        function [LHS] = assembleLHS(obj,s)
            bcapp = obj.bcApplier;
            CtDir = bcapp.computeLinearConditionsMatrix('Dirac');
            CtPer = bcapp.computeLinearPeriodicConditionsMatrix();
            Ct = [CtPer; CtDir];
            K  = s.stiffness;
            C   = Ct';
            nC  = size(Ct,1);
            Z   = zeros(nC);
            LHS = [K C; C' Z];
        end
        
        function [RHS] = assembleRHS(obj,s)
            F = s.forces;
            bcs = obj.boundaryConditions;
            nPer = length(bcs.periodic_leader);
            RHS = [F; zeros(nPer,1); bcs.dirichlet_vals];
        end
        
        function [u,L] = cleanSolution(obj,sol,s)
            nDofsU = size(s.stiffness,1);
            u = sol(1:nDofsU, :);
            L = -sol( (nDofsU+1):end, : );
        end
    end
    
end