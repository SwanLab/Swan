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
    
    methods (Access = public)
        
        function obj = PCG(cParams)
            obj.init(cParams)
            
        end
        
        function x = solve(obj,A,B,x0)
            x = x0;
            r = B - A * x;
            %             z = ModalTesting.applyPreconditioner(M,r);
            z = obj.Preconditioner.apply(r);
            %             z=r-z;
            p = z;
            rzold = r' * z;
            iter = 0;

            Converged = false;

            while not(Converged)
                Ap = A * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                %                 z = ModalTesting.applyPreconditioner(M,r);
                z = obj.Preconditioner.apply(r);
                rznew = r' * z;

                %hasNotConverged = sqrt(rsnew) > tol;
                Converged = norm(r) < obj.tol;

                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
%                 err(iter)=norm(x-xsol);
%                 errAnorm(iter)=((x-xsol)')*A*(x-xsol);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.Preconditioner = cParams.preconditioner;
            obj.tol = cParams.tol;
        end
        
    end
    
end