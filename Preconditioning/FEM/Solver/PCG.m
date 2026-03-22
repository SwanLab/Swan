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
               
        function [x,residual,err,errAnorm] = solve(A,B,x0,P,tol,xsol,mesh,bcApplier,isPlot)
            
            % if nargin == 5, xsol = zeros(size(B)); end
            
            % CANVI: S'HA CANVIAT LA LINIA DE DALT PER AIXÒ
            if nargin < 6, xsol = zeros(size(B)); end
            if nargin < 7, mesh = []; end
            if nargin < 8, bcApplier = []; end
            if nargin < 9, isPlot = false; end
            % FI
            
            iter = 0;
            normB = norm(B);
            x = x0;
%             EIFEMtesting.plotSolution(x,mesh,25,5,iter,bcApplier,0)
            r = B - A(x);
            z = P(r);
            p = z;
            rzold = r' * z;
            normB = norm(B);
            while norm(r)/normB > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                z = P(r);
                rznew = r' * z;
                beta  = (rznew / rzold);
                p = z + beta * p;
                rzold = rznew;

                if isPlot==true
                    x1=bcApplier.reducedToFullVectorDirichlet(x);
                    s.mesh = mesh;
                    s.ndimf = mesh.ndim;
                    s.order = 'P1';
                    s.fValues = reshape(x1,2,[])';
                    name= ['iter' num2str(iter)];
                    fun1=LagrangianFunction(s);
                    fun1.print(name);
                end
                iter = iter + 1;
%               EIFEMtesting.plotSolution(x,mesh,25,5,iter,bcApplier,0)
                residual(iter) = norm(r)/normB;
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A(x-xsol);
            end
        end
        
    end
    
    
end
