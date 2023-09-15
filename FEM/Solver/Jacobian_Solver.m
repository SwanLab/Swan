classdef Jacobian_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS)
            normVal = Inf;
            tol = 1e-6;
            n = length(LHS);
            x = zeros(n,1);
            while normVal>tol
            xold=x;
                for i=1:n
                    sigma=0;
                    for j=1:n
                        if j~=i
                            sigma=sigma+LHS(i,j)*x(j);
                        end
                    end
                    x(i)=(1/LHS(i,i))*(RHS(i)-sigma);
                end
                normVal=norm(xold-x);
            end
        
        end
    end

end