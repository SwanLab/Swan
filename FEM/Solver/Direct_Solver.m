classdef Direct_Solver < Solver
    %UNTITLED25 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        % Analytical Solver (AX = b)
        function x = solve(LHS,RHS)
            x = LHS\RHS;
        end
    end
end

