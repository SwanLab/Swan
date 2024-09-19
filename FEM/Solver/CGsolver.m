classdef CGsolver < handle

    % Cuellos de botella:
    % - [CGSolver]: Line 22. I.e Ku=f.
    % - [Projections and PDE Filters] Systems Ax=b with just 1 dof/node may take few time but nothing compared to the elastic problem matrix system.
    % - [ElasticProblem] Elastic stiffness matrix computation
    % - [ElasticProblem] computeStrain and computeStress takes a little bit
    % - [ComplianceFromConstitutive] Line 47 takes few time, i.e, Integrator.compute

    properties (Access = private)
      x0  
    end

    methods (Access = public)

        function obj = CGsolver()
            obj.init()
        end

        function x = solve(obj,A,b)
            obj.prepareProblem(A);
            tol = 1e-5;
            maxit = 15000;
            tic;
            %P = obj.computePreConditionerMatrix(A,'JACOBI');
            P = obj.computePreConditionerMatrix(A,'ILU0');
            Af = @(x) A*x;            
            x = obj.computeConjugateGradient(Af,b,P,tol,obj.x0,maxit);
            toc;
            obj.x0 = x;
        end

    end

    methods (Access = private)

        function init(obj)
        end

        function prepareProblem(obj, A)
            n = size(A,1);
            if isempty(obj.x0)
                obj.x0 = zeros(n, 1);
            end
        end

        function P = computePreConditionerMatrix(obj,A,type)
            switch type
                case 'IDENTITY'
                    P = @(r) r;
                case 'ILU0'
                    Lchol = ichol(A);
                    P = @(r) obj.applyILU(r,Lchol);
                case 'GAUSS-SEIDEL'
                    P = tril(A);
                case 'JACOBI'
                    Ad = diag(A);
                    P = @(r) r./Ad;
            end
        end        

         function x = computeConjugateGradient(obj,A,b,P,tol,x0,maxit)
            %x = pcg(A,b,tol,maxit,M,M',obj.x0); 
            x = pcg(A,b,tol,maxit,P,[],x0); 
           % x = obj.preconditionedConjugateGradient(A,b,P,tol,x0);
        end


        

    end

    methods (Static, Access = private)


        function [x,residual] = preconditionedConjugateGradient(A,B,P,tol,x0)
            iter = 0;
            n = length(B);
            x = x0;
            r = B - A(x);
            z = P(r);
            p = z;
            rzold = r' * z;
            while norm(r) > tol
                Ap = A(p);
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                z = P(r);
                rznew = r' * z;
                beta  = (rznew / rzold);
                p = z + beta * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r);
            end
        end

       function z = applyILU(r,Lchol)
            z = Lchol\r;
            z = (Lchol')\z;
        end        


    end
end