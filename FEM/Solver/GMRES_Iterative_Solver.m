classdef GMRES_Iterative_Solver < Solver

    methods (Static)
        %GMRES+Arnoldi Iterative Method (Ax = b) (LHSx=RHS)
        function x_k = solve(LHS,RHS)
            k = 5000;
            V = zeros(size(RHS,1),k+1);
            V(:,1) = RHS/norm(RHS);
            H_k = zeros(k+1,k);

            for j = 1:k
                sum = 0;
                z = LHS*V(:,j);
                for i = 1:j
                    H_k(i,j) = z'*V(:,i);
                    sum = sum + H_k(i,j)*V(:,i);
                end
                v = z-sum;
                H_k(j+1,j) = norm(v);
                V(:,j+1) = v/H_k(j+1,j);
            end

            e1 = zeros(k+1,1);
            e1(1,1)=1;

            C = H_k'*H_k;
            D = H_k'*norm(RHS)*e1;
            y_k = C\D;
            x_k = V(:,1:k)*y_k;
        end
    end
end