classdef PCG < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        Preconditioner
        tol
        x0
    end
    
    methods (Static, Access = public)
               
        function [x,residual,err,errAnorm] = solve(A,B,x0,P,tol,xsol,mesh,bcApplier)
            if nargin == 5, xsol = zeros(size(B)); end
            iter = 0;
            normB = norm(B);
            x = x0;
            r = B - A(x);
            z = P(r);
            p = z;
            rzold = r' * z;
            while norm(r)/normB > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                %EIFEMtesting.plotSolution(x,mesh,15,1,iter,bcApplier,0)
                r = r - alpha * Ap;
                z = P(r);
                rznew = r' * z;
                beta  = (rznew / rzold);
                p = z + beta * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r)/normB;
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);
            end
        end
        
    end
    
    
end