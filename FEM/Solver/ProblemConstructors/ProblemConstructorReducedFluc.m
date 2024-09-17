classdef ProblemConstructorReducedFluc < ProblemSolver
    
    properties (Access = private)
        lead, fllw, drch
    end
    
    methods (Access = public)
        
        function obj = ProblemConstructorReducedFluc(cParams)
            obj@ProblemSolver(cParams);
            
            bcs = obj.boundaryConditions;
            obj.lead = bcs.periodic_leader;
            obj.fllw = bcs.periodic_follower;
            obj.drch = bcs.dirichlet_dofs;
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
            [free, ~] = obj.obtainFreeDofs(s.stiffness);
            K = s.stiffness;
            K_II = K(free,free);
            K_IP = K(free,obj.lead) + K(free,obj.fllw); %Grouping P and Q nodal values
            K_PI = K(obj.lead,free) + K(obj.fllw,free); % Adding P  and Q equation
            K_PP = obj.computeStiffnessPeriodicPeriodic(K); % Adding and grouping
            LHS = [K_II, K_IP; K_PI, K_PP];
        end
        
        function [RHS] = assembleRHS(obj,s)
            [free, ~] = obj.obtainFreeDofs(s.stiffness);
            F = s.forces;
            F_I = F(free);
            F_P = F(obj.lead) + F(obj.fllw);
            RHS = [F_I; F_P];
        end
        
        function [u,L] = cleanSolution(obj,sol,s)
            bcs   = obj.boundaryConditions;
            [free, nDofs] = obj.obtainFreeDofs(s.stiffness);
            u = zeros(nDofs,1);
            u(free) = sol(1:1:size(free,2));
            u(obj.lead) = sol(size(free,2)+1:1:size(sol,1));
            u(obj.fllw) = u(obj.lead);
            u(obj.drch) = bcs.dirichlet_vals;
            L = [];
        end
        
        function [free, nDofs] = obtainFreeDofs(obj,K)
            nDofs = size(K,1);
            dofs = 1:nDofs;
            free = setdiff(dofs, [obj.lead; obj.fllw; obj.drch]);
        end
        
        function [K_PP] = computeStiffnessPeriodicPeriodic(obj,K)
            K_LL = K(obj.lead,obj.lead); 
            K_LF = K(obj.lead,obj.fllw);
            K_FL = K(obj.fllw,obj.lead);
            K_FF = K(obj.fllw,obj.fllw);
            K_PP = K_LL + K_LF + K_FL + K_FF;
        end
    end
    
end