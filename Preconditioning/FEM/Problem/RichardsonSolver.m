classdef RichardsonSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Static, Access = public)
          

        function [x,residual,err,errAnorm] = solve(A,B,x0,P,tol,linesearch,xsol,mesh,bcApplier)
            if nargin == 6, xsol = zeros(size(B)); end            
            iter = 0;
            n = length(B);
            x = x0;
            r = B - A(x);      
            z = P(r);
            while norm(r) > tol
%                 tau = 1*linesearch(r,A);
                 tau = 1*linesearch(z,r,A);
                x = x + tau * z;
%                EIFEMtesting.plotSolution(x,mesh,10,10,iter,bcApplier,0)
                r = B - A(x); 
                z = P(r);
                test(iter+1)=norm(z);
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);  
                %norm(r)
            end        

        end

    end

    
end