classdef PCG < handle

    properties (Access = private)
        Preconditioner
        tol
    end

    properties (Access = private)
        x0
        iter0
    end

    methods (Access = public)
        function obj = PCG(cParams)
            obj.init(cParams);
        end

        function [x,residual,err,errAnorm] = solve(obj,LHS,B,xsol)
            if nargin == 3, xsol = zeros(size(B)); end
            obj.prepareProblem(B);
            P = @(A,r) obj.Preconditioner.solve(A,r);
            A = @(x) LHS*x;
            iter = 0;
            x = obj.x0;
            r = B - A(x);
            z = P(LHS,r);
            p = z;
            rzold = r' * z;
            normB = norm(B);
            while norm(r)/normB > obj.tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                z = P(LHS,r);
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
            obj.checkStep(iter);
            disp(['Iter: ',num2str(iter)]);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.Preconditioner = cParams.preconditioner;
            obj.tol            = cParams.tol;
        end

        function prepareProblem(obj, b)
            n = length(b);
            if isempty(obj.x0)
                obj.x0 = zeros(n, 1);
            end
        end

        function checkStep(obj,iter)
            if isempty(obj.iter0)
                obj.iter0 = iter;
            elseif iter>3*obj.iter0
                obj.iter0 = [];
                obj.Preconditioner.update();
            end
        end
    end
end