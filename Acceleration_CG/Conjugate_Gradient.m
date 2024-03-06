classdef Conjugate_Gradient < handle

    properties
        xOld
        rhoOld
        
        tolMax
        tolMin
        tolStandard
        tol
        rhoNormOld

        globalIters

        designVariable,volume

        rhoNormFactor, tolNormFactor
    end

    methods

        function obj = Conjugate_Gradient(cParams)
            obj.designVariable = cParams.designVariable;
            obj.volume         = cParams.volume;
            obj.setToleranceFunctions();

            obj.globalIters = 1; 
            obj.rhoNormFactor = 4/log10(numel(obj.designVariable.fun.fValues));            
        end

    end
    
    methods (Access = public)

        function x = solve(obj,A,b)
            tic
            % Design variable
            rho        = obj.designVariable.fun.fValues;
            rhoNorm    = obj.returnDesignVarIncNorm(rho);
            obj.rhoOld = rho;
            if isempty(obj.xOld)
                x = A\b;                
                disp('Iter: DIRECT')
                obj.tolNormFactor = 4.5/norm(b);
                obj.tol = 1e-1*obj.tolNormFactor;
                obj.rhoNormOld = 0;
            else
                obj.obtainTolerance(rhoNorm);
                x   = obj.xOld;
                [x,~,~,it] = pcg(A,b,obj.tol,3e3,[],[],x);
                disp('Iter: ' + string(it))
                disp('Tol: ' + string(obj.tol/obj.tolNormFactor))
            end
            % disp('E: ' + string(1/2*x'*A*x - b'*x))
            disp('rho normInc: ' + string(rhoNorm))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
            obj.xOld = x;
        end

        function rhoNorm = returnDesignVarIncNorm(obj,rho)
            if isempty(obj.rhoOld)
                rhoNorm = 0;
            else
                rhoNorm = norm(rho - obj.rhoOld);
            end
        end

    end

    methods (Access = private)

        function obtainTolerance(obj,rhoNorm)
            d = obj.designVariable;
            [c,~] = obj.volume.computeFunctionAndGradient(d);
            if rhoNorm >= obj.rhoNormOld
                if rhoNorm > 1/obj.rhoNormFactor
                    obj.tol = obj.tolStandard(c)*obj.rhoNormFactor*rhoNorm;
                elseif rhoNorm > 0 && rhoNorm < 1
                    obj.tol = 0.1*obj.tolStandard(c);
                end
                obj.tol = min(obj.tolMax(c),max(obj.tolMin(c),obj.tol))*obj.tolNormFactor;
            end
            obj.rhoNormOld = rhoNorm;
        end

        function setToleranceFunctions(obj)
            d = obj.designVariable;
            [cInit,~] = obj.volume.computeFunctionAndGradient(d);
            obj.tolStandard = @(c) interp1([0,0.2,1],[5e-3,1e-2,7e-2],abs(c/cInit));
            obj.tolMax      = @(c) interp1([0,0.8,1],[7e-2,1e-1,5e-2],abs(c/cInit));          
            obj.tolMin      = @(c) interp1([0,0.5,1],[1e-4,1e-3,3e-3],abs(c/cInit));
            % obj.tolStandard = @(c) interp1([0,0.2,1],[1e-3,5e-3,1e-3],abs(c/cInit));
            % obj.tolMax      = @(c) interp1([0,0.8,1],[1e-2,5e-2,1e-2],abs(c/cInit));          
            % obj.tolMin      = @(c) interp1([0,0.5,1],[1e-5,5e-4,1e-4],abs(c/cInit));
        end

    end

    methods (Static, Access = private)

    end
end