classdef ConjugateGradient < handle

    properties (Access = private)
        xOld
        tol
        maxIters
        displayInfo
    end

    methods (Access = public)

        function obj = ConjugateGradient(cParams)
            obj.init(cParams);
        end

        function x = solve(obj,A,b)
            tic
            if isempty(obj.xOld) % - Might be replaced in the future by an stored initial guess
                x  = A\b;
                it = -1;
            else
                x = obj.xOld;
                t = obj.tol.val;                
                [x,~,~,it] = pcg(A,b,t,obj.maxIters,[],[],x);
            end
            obj.xOld = x;            
            obj.displayInfo(it);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.tol      = cParams.tol;
            obj.maxIters = cParams.solverParams.maxIters;
            if cParams.solverParams.displayInfo
                obj.displayInfo = @obj.printSolverInfo;
            else
                obj.displayInfo = @obj.emptyFunc;
            end
        end
        
        function printSolverInfo(obj,it)
            disp('Tol: ' + string(obj.tol.val))
            disp('Iter: ' + string(it))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
        end

    end

    methods (Static, Access = private)
        
        function emptyFunc(~)
            
        end

    end

end