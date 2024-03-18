classdef ConjugateGradient < handle

    properties (Access = private)
        xOld
        tol
        maxIters
        tolFactor
        xBest
    end

    methods (Access = public)

        function obj = ConjugateGradient(cParams)
            obj.tol      = cParams.tol;
            obj.maxIters = 3e3; % - 
        end

        function x = solve(obj,A,b)
            tic
            if isempty(obj.xOld) % - Might be replaced in the future by an stored initial guess
                obj.defineCorrectionFactor(b);
                x         = A\b;
                obj.xBest = x;
                disp('DIRECT METHOD')
            else
                t = obj.tol.val*obj.tolFactor;
                rBest = norm(b - A*obj.xBest);
                rOld  = norm(b - A*obj.xOld);
                if rBest <= rOld
                    x = obj.xBest;
                else
                    x = obj.xOld;
                    obj.xBest = obj.xOld;
                end
                [x,~,~,it] = pcg(A,b,t,obj.maxIters,[],[],x);
                disp('Iter: ' + string(it))
            end
            obj.xOld = x;
            disp('Tol: ' + string(obj.tol.val))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
        end

    end

    methods (Access = private)

        function defineCorrectionFactor(obj,b)
            obj.tolFactor = 4/norm(b);
        end

    end
end