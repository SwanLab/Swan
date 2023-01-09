classdef OptimizerNullSpace < Optimizer

    properties (Access = public)
        aGtNum = 1.5e-3
        aJ  = 1
        DeltaNum = 100000
        Delta
        aGt
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
            obj.aGt   = obj.aGtNum/(abs(obj.constraint.value)+0.1)^1;
%             obj.Delta = obj.DeltaNum/(abs(obj.constraint.value)+sqrt(0.01))^2;
%             obj.aGt = obj.aGtNum;
            obj.Delta = obj.DeltaNum;

            x0 = obj.designVariable.value;
            obj.saveOldValues(x0);
            if obj.nIter == 0
                obj.calculateInitialStep(); % for merit old (?)
            end
            obj.dualUpdater.t = obj.primalUpdater.tau;
            obj.dualUpdater.update();
            s.x = x0;
            tOld = obj.primalUpdater.tau;
            s.t = tOld;
            obj.mOld = obj.computePreviousMeritFunction(s);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;
            while ~obj.acceptableStep
                obj.computeMeritGradient(DJ,Dg,tOld);
                x = obj.updatePrimal();
                obj.designVariable.update(x);
                obj.checkStep(x,x0,tOld);
            end
            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            if obj.nIter == 0
                obj.primalUpdater.computeFirstStepLength(1);
            else
                factor = 1.2;
                obj.primalUpdater.increaseStepLength(1);
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

        function computeMeritGradient(obj,DJ,Dg,t)
            lJ  = obj.dualUpdater.lJ;
            lG  = obj.projectLambdaG();
            f   = obj.aGt/(t*obj.aJ);
%             l   = lJ+f*lG;
            l   = t*(lJ+f*lG);
            DmF = obj.aJ*(DJ + Dg*l);
            obj.meritGradient = DmF;
        end

        function checkStep(obj,x,x0,t0)
            s.x = x;
            s.t = obj.primalUpdater.tau;
            mNew = obj.computeMeritFunction(s);
            if mNew < obj.mOld && norm(x-x0)/norm(x0) < 10
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            elseif obj.primalUpdater.isTooSmall()
%                 error('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew = obj.mOld;
                obj.primalUpdater.increaseStepLength(t0); % restore tau to prev step
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function mF = computeMeritFunction(obj,cParams)
            x  = cParams.x;
            t  = cParams.t;
            lJ = obj.dualUpdater.lJ;
            lG = obj.projectLambdaG();
            obj.designVariable.update(x)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J  = obj.cost.value;
            h  = obj.constraint.value;
            aG = obj.aGt/t;
%             mJ = obj.aJ*(J+lJ*h);
            mJ = obj.aJ*(J+t*lJ*h);
%             mG = aG*lG*h;
            mG = aG*t*lG*h;
            mF = mJ+mG;
        end

        function mF = computePreviousMeritFunction(obj,cParams)
            mF = obj.computeMeritFunction(cParams);
            if obj.nIter == 0
                mF = 1*mF; % "Initial Lagrange multiplier"
            end
        end

        function lG = projectLambdaG(obj)
            lG = obj.dualUpdater.lG;
            lG = max(-obj.Delta,min(obj.Delta,lG));
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