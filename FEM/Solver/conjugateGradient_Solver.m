classdef conjugateGradient_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS)
            tol = 1e6;
            n = length(RHS);
            x = zeros(n, 1); 
            r = RHS - LHS * x; 
            p = r; 
            rsold = r' * r;

            hasNotConverged = true;

            while hasNotConverged
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;

                hasNotConverged = sqrt(rsnew) < tol;
                

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            end
        
        end
    end
end