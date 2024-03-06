classdef Direct_Solver < Solver

    methods (Static)
        % Analytical Solver (AX = b)
        
        function x = solve(LHS,RHS)
            tic
            x = LHS\RHS;
            disp('DIRECT SOLVER')
            disp('Convergence time ' + string(toc) + ' s')
        end
 
    end

end