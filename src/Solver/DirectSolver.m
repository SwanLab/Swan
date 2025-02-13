classdef DirectSolver < Solver

    methods (Static)
        function x = solve(LHS,RHS)
            x = LHS\RHS;
        end
    end

end