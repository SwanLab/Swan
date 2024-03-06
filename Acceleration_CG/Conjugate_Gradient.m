classdef Conjugate_Gradient < handle

    properties
        xOld
        rhoOld
        
        tolMax
        tolMin
        tolStandard

        globalIters

        designVariable,volume

        rhoNormFactor, tolNormFactor
    end

    methods

        function obj = Conjugate_Gradient(cParams)
            obj.designVariable = cParams.designVariable;
            obj.volume         = cParams.volume;
            obj.setToleranceFunction();

            obj.globalIters = 1; 
            obj.tolMax      = 1e-1;
            obj.tolMin      = 1e-5;

            obj.rhoNormFactor = 4/log10(numel(obj.designVariable.fun.fValues));
            
        end

    end
    
    methods (Access = public)

        function x = solve(obj,A,b)
            tic
            % Design variable
            rho        = obj.designVariable.fun.fValues;
            rhoNorm    = obj.returnDesignVarNorm(rho);
            obj.rhoOld = rho;
            if isempty(obj.xOld)
                x = A\b;                
                disp('Iter: DIRECT')
                obj.tolNormFactor = 4.5/norm(b);
            else
                tol = obj.obtainTolerance(rhoNorm);
                x   = obj.xOld;
                [x,~,~,it] = pcg(A,b,tol,3e3,[],[],x);
                disp('Iter: ' + string(it))
                disp('Tol: ' + string(tol/obj.tolNormFactor))
            end
            disp('E: ' + string(1/2*x'*A*x - b'*x))
            disp('rho normInc: ' + string(rhoNorm))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
            obj.xOld = x;
        end

        function rhoNorm = returnDesignVarNorm(obj,rho)
            if isempty(obj.rhoOld)
                rhoNorm = 0;
            else
                rhoNorm = norm(rho - obj.rhoOld);
            end
        end

    end

    methods (Access = private)

        function tol = obtainTolerance(obj,rhoNorm)
            d = obj.designVariable;
            [c,~] = obj.volume.computeFunctionAndGradient(d);
            if rhoNorm > 1/obj.rhoNormFactor
                tol = obj.tolStandard(c)*obj.rhoNormFactor*rhoNorm;%(1 + rhoNorm/2);
            else
                tol = 0.1*obj.tolStandard(c);
            end
            tol = min(obj.tolMax,max(obj.tolMin,tol))*obj.tolNormFactor;
        end

        function setToleranceFunction(obj)
            d = obj.designVariable;
            [cInit,~] = obj.volume.computeFunctionAndGradient(d);
            obj.tolStandard = @(c) interp1([0,0.2,1],[1e-3,1e-2,7e-2],abs(c/cInit));
        end

    end
end