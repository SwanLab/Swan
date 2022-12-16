classdef OptimizerNullSpace < Optimizer

    properties (Access = public)
        aGt = 1/sqrt(10)  % 1e-2
        aJ  = sqrt(10)    % 1e2
    end

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
        fProv
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        costOld
        upperBound
        lowerBound
        tol = 1e-8
        nX
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        incrementalScheme
        hasFinished
        mOld
        meritNew
        nConstr
        meritGradient

        globalCost
        globalConstraint
        globalCostGradient
        globalMerit
        globalLineSearch
        globalDual
        globalDesignVar
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams,optParams)
            obj.fProv = optParams.fProv;
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);

            obj.dualUpdater.parameter = optParams.fProv;

            obj.prepareFirstIter();
        end

        function obj = solveProblem(obj)
            obj.hasConverged = false;
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.hasFinished = false;
            obj.printOptimizerVariable();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.printOptimizerVariable();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.upperBound             = cParams.uncOptimizerSettings.ub;
            obj.lowerBound             = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.designVariable         = cParams.designVar;
            obj.dualVariable           = cParams.dualVariable;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.nConstr                = cParams.constraint.nSF;
            obj.nX                     = length(obj.designVariable.value);
            obj.maxIter                = cParams.maxIter;
            obj.hasConverged           = false;
            obj.nIter                  = 0;
        end

        function prepareFirstIter(obj)
            obj.cost.computeFunctionAndGradient();
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.value = zeros(obj.nConstr,1);
        end

        function obj = update(obj)
            x0 = obj.designVariable.value;
            obj.saveOldValues(x0);
            if obj.nIter ~= 0
                obj.dualUpdater.update(); 
            end
            obj.calculateInitialStep();
            obj.mOld = obj.computeMeritFunction(x0);
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;
            while ~obj.acceptableStep

                obj.dualUpdater.parameter = obj.aGt/(obj.aJ*obj.primalUpdater.tau); % NEW NEW NEW for EQUALITY

                obj.computeMeritGradient(DJ,Dg);
                x = obj.updatePrimal();
                obj.designVariable.update(x);
                obj.dualUpdater.update();
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            if obj.nIter == 0
                obj.cost.computeFunctionAndGradient();
                obj.constraint.computeFunctionAndGradient();
                x       = obj.designVariable.value;
                l       = obj.dualVariable.value;
                DJ      = obj.cost.gradient;
                Dg      = obj.constraint.gradient;
                DmF     = obj.aJ*(DJ + Dg*l);
                factor  = 1;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor,1);
            else
                factor = 1.2; % 1.2
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function displayIter(obj,x)
            m = obj.designVariable.mesh.innerMeshOLD;
            bm = m.createBoundaryMesh();
            s.backgroundMesh = m;
            s.boundaryMesh   = bm;
            um = UnfittedMesh(s);
            um.compute(x);
            figure()
            um.plot();
        end

        function x = updatePrimal(obj)
            x       = obj.designVariable.value;
            g       = obj.meritGradient;
            x       = obj.primalUpdater.update(g,x);
        end

        function computeMeritGradient(obj,DJ,Dg)
            l       = obj.dualVariable.value;
            DmF     = obj.aJ*(DJ + Dg*l);
            obj.meritGradient = DmF;
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);
            if obj.nIter == 0 %&& mNew <= obj.mOld
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
                obj.dualUpdater.updateOld();
            end
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            elseif obj.primalUpdater.isTooSmall()
%                 error('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function mF = computeMeritFunction(obj,x)
            DJ = obj.cost.gradient;
            Dh = obj.constraint.gradient;
            S  = (Dh'*Dh)^-1;
            lJ = -S*(Dh'*DJ);
            obj.designVariable.update(x)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J  = obj.cost.value;
            h  = obj.constraint.value;
            t = obj.primalUpdater.tau;
            aG = obj.aGt/t;
            mG = 0.5*aG*h'*S*h;
            mJ = obj.aJ*(J + lJ'*h);
            mF = mJ+mG;
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
        end

        function obj = updateOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
        end

        function obj = checkConvergence(obj)
           if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
               obj.hasConverged = true;
           else
               
           end

        end

        function obj = updateMonitoring(obj)
            s.nIter            = obj.nIter;
            s.tau              = obj.primalUpdater.tau;
            s.lineSearch       = obj.lineSearch;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.oldCost          = obj.oldCost;
            s.hasFinished      = obj.hasFinished;
            s.meritNew         = obj.meritNew;
            obj.outputFunction.monitoring.compute(s);
        end

        function updateIterInfo(obj)
            obj.increaseIter();
            obj.updateStatus();
        end

        function increaseIter(obj)
            obj.nIter = obj.nIter + 1;
        end

        function updateStatus(obj)
            obj.hasFinished = obj.hasConverged || obj.hasExceededStepIterations();
        end

        function itHas = hasExceededStepIterations(obj)
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end

        function saveVariablesForAnalysis(obj)
            i                           = obj.nIter + 1;
            obj.globalCost(i)           = obj.cost.value;
            obj.globalConstraint(:,i)   = obj.constraint.value;
            obj.globalCostGradient(i)   = norm(obj.cost.gradient);
%             obj.globalMerit(i)          = obj.meritNew;
%             obj.globalLineSearch(i)     = obj.primalUpdater.tau;
            obj.globalDual(:,i)         = obj.dualVariable.value;
            obj.globalDesignVar(:,i)    = obj.designVariable.value;
            if obj.hasConverged
                c = obj.globalCost;
                h = obj.globalConstraint;
                g = obj.globalCostGradient;
%                 m = obj.globalMerit;
%                 t = obj.globalLineSearch;
                d = obj.globalDual;
                v = obj.globalDesignVar;
                save('NullSpaceCant04.mat',"c","g","h","d","v");
            end
        end

    end

end