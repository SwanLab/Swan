classdef Solver
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function d_u = analytical(d_u,LHS,RHS,vR,vL,displacements)
            if ~isempty(vR)
                d_u(vR) = displacements(:,3);
                d_u(vL,1) = LHS(vL,vL)\(RHS(vL) - LHS(vL,vR)*d_u(vR));
            else
                d_u(vL,1) = LHS(vL,vL)\RHS(vL);
            end
        end
    end
    
end

