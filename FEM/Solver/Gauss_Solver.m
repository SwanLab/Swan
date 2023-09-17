classdef Gauss_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS)
            normVal = Inf;
            tol = 1e-6;
            n = length(LHS);
            x = zeros(n,1);
            numItr = 0;
            while normVal>tol
                x_old=x;
                
                for i=1:n
                    
                    sigma=0;
                    
                    for j=1:i-1
                            sigma=sigma+LHS(i,j)*x(j);
                    end
                    
                    for j=i+1:n
                            sigma=sigma+LHS(i,j)*x_old(j);
                    end
                    
                    x(i)=(1/LHS(i,i))*(RHS(i)-sigma);
                end
                
                numItr=numItr+1;
                normVal=norm((x_old-x)/x);
            end
        end
    end

end