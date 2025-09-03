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

    methods (Access = private)
        function prepareProblem(obj, b)
            n = length(b);
            if isempty(obj.x0)
                obj.x0 = zeros(n, 1);
            end
        end
    end

    methods (Access = public)

        function [x,residual,err,errAnorm] = solve(obj,A,B,P,tol,xsol)
            if nargin == 5, xsol = zeros(size(B)); end
            obj.prepareProblem(B);
            iter = 0;
            x = obj.x0;
            r = B - A(x);
            z = P(r);
            p = z;
            rzold = r' * z;
            normB = norm(B);
            while norm(r)/normB > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
%                  EIFEMtesting_3D.plotSolution(x,mesh,30,4,iter,bcApplier,0)
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
            obj.x0 = x;
            disp(['Iter: ',num2str(iter)]);
        end
        
    end
    
    
end