classdef Solver_Analytical<Solver
    %AnalyticalSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = public,Static)
        
        % Analytical Solver (A·X=b)
        function x = solve(LHS,RHS,dof,fixnodes)
            x = zeros(dof.ndof,1);
            if ~isempty(dof.vR)
                x(dof.vR) = fixnodes(:,3);
                x(dof.vL,1) = LHS(dof.vL,dof.vL)\(RHS(dof.vL) - LHS(dof.vL,dof.vR)*x(dof.vR));
            else
                x(dof.vL,1) = LHS(dof.vL,dof.vL)\RHS(dof.vL);
            end
        end
    end
    
end

