classdef Direct_Solver < Solver

    methods (Static)
        % Analytical Solver (AX = b)
        function x = solve(LHS,RHS)
            x = LHS\RHS;
        end
    end

end