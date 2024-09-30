classdef RichardsonSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = RichardsonSolver()
            
        end

        function [x,residual,err,errAnorm] = solve(A,B,x0,P,tol,tau,xsol)
            if nargin == 6, xsol = zeros(size(B)); end            
            iter = 0;
            n = length(B);
            x = x0;
            r = B - A(x);      
            z = P(r);
            while norm(r) > tol
                x = x + tau * z;
                r = B - A(x); 
                z = P(r);
                test(iter+1)=norm(z);
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);                
            end        

        end
    end

    
end