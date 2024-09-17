classdef ProblemConstructorMonolithicDisp < ProblemSolver

    properties (Access = private)
        bcApplier
        notPeriodic
    end
    
    methods (Access = public)

        function obj = ProblemConstructorMonolithicDisp(cParams)
            obj@ProblemSolver(cParams);
            obj.bcApplier = cParams.BCApplier;
            
            bc = obj.boundaryConditions;
            obj.notPeriodic = isempty(bc.periodic_leader);
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
            if obj.notPeriodic
                Ct = bcapp.computeLinearConditionsMatrix('Dirac');
            else
                iV = s.iVoigt;
                nV = s.nVoigt;
                Ct = bcapp.computeSingleDirichletPeriodicCondition(iV, nV);
            end
            K   = s.stiffness;
            C   = Ct';
            nC  = size(Ct,1);
            Z   = zeros(nC);
            LHS = [K C; C' Z];
        end
        
        function [RHS] = assembleRHS(obj,s)
            bcs = obj.boundaryConditions;
            bcapp = obj.bcApplier;
            if obj.notPeriodic
                F      = s.forces;
                nCases = size(F,2);
                lambda = bcs.dirichlet_vals;
                Ct     = repmat(lambda, [1 nCases]);
                RHS = [F; Ct];
            else
                iV = s.iVoigt;
                nV = s.nVoigt;
                RHS = bcapp.computeMicroDisplMonolithicRHS(iV, nV);
            end
        end
        
        function [u,L] = cleanSolution(obj,sol,s)
            nDofsU = size(s.stiffness,1);
            u = sol(1:nDofsU, :);
            L = -sol( (nDofsU+1):end, : );
        end
    end
    
end