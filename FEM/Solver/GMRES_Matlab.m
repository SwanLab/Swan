classdef GMRES_Matlab

    methods (Static)

        function x_k = solve(LHS,RHS)
            A = LHS;
            b = RHS;
            MAXIT = 5000;
%             x_k = gmres(A,b,[],[],MAXIT);
            x_k = minres(A,b,[],MAXIT);
        end
    end
end

