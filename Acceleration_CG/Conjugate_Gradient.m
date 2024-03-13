classdef Conjugate_Gradient < handle

    properties (Access = private)
        xOld, xStep
        rhoIter, rhoOld, rhoOld_step
        rhoNormIter, rhoNormOld, rhoNorm_step
        tolMax, tolMin, tolStandard, tol
        rhoNormFactor, tolNormFactor
        isZero
        designVariable,volume
        defineToleranceValue
    end

    methods (Access = public)

        function obj = Conjugate_Gradient(cParams)
            obj.init(cParams);          
        end

        function x = solve(obj,A,b)
            tic
            obj.updateDesignVarIter();
            if isempty(obj.xOld)
                x = A\b;                
                disp('Iter: DIRECT')
                obj.computeFirstIterParameters(b);
                obj.xOld = x;
            else
                obj.defineToleranceValue();
                x = obj.xOld;
                [x,~,~,it] = pcg(A,b,obj.tol,3e3,[],[],x);
                disp('Iter: ' + string(it))
                disp('Tol: '  + string(obj.tol/obj.tolNormFactor))
            end
            disp('rho normInc (scaled): ' + string(obj.rhoNormIter*obj.rhoNormFactor))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
            % obj.xAccepted = obj.xOld;
            obj.xOld = x;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.solverVars.designVariable;
            obj.volume         = cParams.solverVars.volume;
            obj.setToleranceDefinitionStrategy(cParams);            
            obj.rhoNormFactor = 4/log10(numel(obj.designVariable.fun.fValues)); % Scaled with number of design vars  
        end

        function computeProjectedGradientTolerance(obj)
            rNIter = obj.computeDesignVarIncNorm(obj.rhoIter,obj.rhoOld);
            obj.rhoNormIter = rNIter;
            if rNIter > 0
                % - ||rho_k+1 - rho_k|| > 0 -
                rhoNormStep = obj.computeDesignVarIncNorm(obj.rhoIter,obj.rhoOld_step);
                if rhoNormStep > 0
                    if obj.isZero
                        % - Increasing step -
                        disp('- New step... -')
                        obj.updateDesignVarPreviousStep(obj.rhoOld);
                        d       = obj.designVariable;
                        [c,~]   = obj.volume.computeFunctionAndGradient(d);
                        t       = obj.tolStandard(c)*max(1,obj.rhoNormFactor*obj.rhoNormIter);
                        obj.tol = min(obj.tolMax(c),max(obj.tolMin(c),t))*obj.tolNormFactor;
                        obj.isZero = false;
                    else
                        % - Decreasing step -
                        disp('- Decreasing step... -')
                    end
                    obj.updateDesignVarPrevious(obj.rhoIter);
                    obj.updateOldNorm(rhoNormStep);
                else
                    % - Restarting merit function step -
                    disp('- Restarting merit function... ')
                    obj.tol = 1e-3;
                end
            else
                % - Accepted step, recomputing gradient (?) -
                obj.isZero = true;
            end
        end

        function computeMMATolerance(obj)
            rNIter          = obj.computeDesignVarIncNorm(obj.rhoIter,obj.rhoOld);
            obj.rhoNormIter = rNIter;
            rNScaled        = rNIter*obj.rhoNormFactor;
            t               = obj.tolStandard(rNScaled)*rNScaled;
            obj.tol         = min(obj.tolMax,max(obj.tolMin,t))*obj.tolNormFactor;
            obj.updateDesignVarPrevious(obj.rhoIter);
        end

        function setToleranceDefinitionStrategy(obj,cParams)
            switch cParams.solverVars.optimizer
                case 'MMA'
                    obj.defineToleranceValue = @obj.computeMMATolerance;
                    obj.setMMAToleranceBounds();
                case 'Null Space'
                    obj.defineToleranceValue = @obj.computeProjectedGradientTolerance;
                    obj.setProjectedGradientToleranceBounds();
                otherwise
                    error('Tolerance update strategy not implemented for this optimizer')
            end
        end

        function setProjectedGradientToleranceBounds(obj)
            d = obj.designVariable;
            [cInit,~] = obj.volume.computeFunctionAndGradient(d);
            obj.tolStandard = @(c) interp1([0,0.2,1],[1e-3,1e-2,5e-2],abs(c/cInit));
            obj.tolMax      = @(c) interp1([0,0.2,1],[1e-2,6e-2,2e-1],abs(c/cInit));
            obj.tolMin      = @(c) interp1([0,0.5,1],[1e-5,5e-4,5e-4],abs(c/cInit));
        end

        function setMMAToleranceBounds(obj)
            obj.tolStandard = @(rhoNorm) interp1([0,0.01,0.05,0.21,1,1e5],[1e-12,1e-5,1e-5,1e-1,10,10],rhoNorm);
            obj.tolMax      = 8e-1;
            obj.tolMin      = 1e-8;
        end

        function updateDesignVarIter(obj)            
            obj.rhoIter = obj.designVariable.fun.fValues;
        end

        function updateDesignVarPreviousStep(obj,rho)
            obj.rhoOld_step = rho;
        end

        function updateDesignVarPrevious(obj,rho)
            obj.rhoOld = rho;
        end

        function updateOldNorm(obj,rhoNorm)
            obj.rhoNormOld = rhoNorm;
        end

        function relation = computeDesignVarIncNormRelation(obj,rhoNorm)
            relation = rhoNorm/obj.rhoNormOld;
        end

        function computeFirstIterParameters(obj,b)
            obj.tolNormFactor = 4.5/norm(b); % Factor for pcg function
            obj.tol           = 0.1*obj.tolNormFactor;
            obj.updateDesignVarPrevious(obj.rhoIter);
            obj.updateDesignVarPreviousStep(obj.rhoIter);
            obj.rhoNormIter  = 1e-5;
            obj.rhoNorm_step = 1e-5;
            obj.rhoNormOld   = 1e-5;
            obj.isZero       = false;
        end

        function setToleranceFunctions(obj)
            d = obj.designVariable;
            [cInit,~] = obj.volume.computeFunctionAndGradient(d);
            obj.tolStandard = @(c) interp1([0,0.2,1],[1e-3,1e-2,1e-1],abs(c/cInit));
            obj.tolMax      = @(c) interp1([0,0.2,1],[1e-2,6e-2,2e-1],abs(c/cInit));          
            obj.tolMin      = @(c) interp1([0,0.5,1],[1e-5,5e-4,5e-4],abs(c/cInit));
        end

    end

    methods (Static, Access = private)

        function rhoNorm = computeDesignVarIncNorm(rho1,rho2)
            rhoNorm = norm(rho1 - rho2);
        end
        
    end

end