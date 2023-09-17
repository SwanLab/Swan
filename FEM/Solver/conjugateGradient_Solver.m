classdef conjugateGradient_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS)
            maxIter = 1e6;
            tol = 1e6;
            n = length(RHS);
            x = zeros(n, 1); % Initial guess
            r = RHS - LHS * x; % Initial residual
            p = r; % Initial search direction
            rsold = r' * r;

            for k = 1:maxIter
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;
        
                % Check for convergence
                if sqrt(rsnew) < tol
                    fprintf('Conjugate Gradient method converged after %d iterations.\n', k);
                    return;
                end
        
                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
            end
        
            fprintf('Conjugate Gradient method did not converge after %d iterations.\n', maxIter);
        end
    end
end