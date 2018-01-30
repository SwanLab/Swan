classdef Solver_Dirichlet_Conditions < Solver
    %AnalyticalSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fixnodes
    end
    
    methods (Access = public)
        
        function obj = Solver_Dirichlet_Conditions(obj)
        end
        %     methods (Access = public,Static)
        % Analytical Solver (A·X=b)
        function x = solve(obj,x,LHS,RHS,dof)
%             x = zeros(dof.ndof,1);
            if ~isempty(dof.vR)
                x(dof.vR) = obj.fixnodes(:,3);
                x(dof.vL,1) = LHS(dof.vL,dof.vL)\(RHS(dof.vL) - LHS(dof.vL,dof.vR)*x(dof.vR));
            else
                x(dof.vL,1) = LHS(dof.vL,dof.vL)\RHS(dof.vL);
            end
        end
        
        function setSolverVariables(obj,data)
            obj.fixnodes = data.fixnodes;
        end
    end
    
end

