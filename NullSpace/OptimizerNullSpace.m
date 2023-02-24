classdef OptimizerNullSpace < Optimizer
% Why for it>18 lG oscillates around negative values?
% - Let's try Florian's merit
% - ...

    properties (Access = public)
        aG
        eta
        p
        n0 = 100
    end

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';

        lGtrial
        lGmax
        lG
        lJ
        l
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
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
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
            if obj.nIter<=obj.n0 % if -> clean code
                obj.p     = 0;
                obj.aG    = 1;
                if obj.nIter == 0
                    obj.eta = inf;
                else
                    obj.eta = 0.02;
                end
            else
                obj.p     = inf;
                obj.aG    = 1;
                obj.eta   = 0.005;
            end
            x0 = obj.designVariable.value;
            g0 = obj.constraint.value;
            obj.saveOldValues(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;

            while ~obj.acceptableStep
                obj.designVariable.update(x0);
                obj.dualUpdater.t  = obj.primalUpdater.tau;
                obj.dualUpdater.aG = obj.aG;
                obj.dualUpdater.p  = obj.p;
                obj.dualUpdater.update();
                obj.mOld = obj.computeMeritFunction(x0);
%                 obj.mOld = obj.computeMeritFunctionFlorian(x0,x0);
                obj.computeMeritGradient(DJ,Dg);
                x = obj.updatePrimal();
                obj.designVariable.update(x);
                obj.checkStep(x,x0,g0);
            end

            obj.lGtrial(obj.nIter+1) = obj.dualUpdater.lGtrialPl;
            obj.lGmax(obj.nIter+1) = obj.dualUpdater.lGmaxPl;
            obj.lG(obj.nIter+1) = obj.dualUpdater.lGPl;
            obj.lJ(obj.nIter+1) = obj.dualUpdater.lJPl;
            obj.l(obj.nIter+1) = obj.dualUpdater.lPl;

            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            if obj.nIter == 0
                obj.primalUpdater.computeFirstStepLength(1);
            else
                factor = 2;
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
            l   = obj.dualVariable.value;
            DmF = DJ+l*Dg;
            obj.meritGradient = DmF;
        end

        function checkStep(obj,x,x0,g0)
            mNew = obj.computeMeritFunction(x);
            g  = obj.constraint.value;
            v0 = obj.computeVolume(g0);
            v  = obj.computeVolume(g);
            if mNew < obj.mOld && norm(v-v0) < obj.eta
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            elseif obj.primalUpdater.isTooSmall()
%                 error('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew = obj.mOld;
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function v = computeVolume(obj,g)
            targetVolume = obj.targetParameters.Vfrac;
            v            = targetVolume*(1+g);
        end

        function mF = computeMeritFunction(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            l  = obj.dualVariable.value;
            J  = obj.cost.value;
            h  = obj.constraint.value;
            mF = J+l*h;
        end

%         function mF = computeMeritFunctionFlorian(obj,x,y)
%             mF = obj.computeMeritFunction(x);
% 
%             obj.designVariable.update(y);
%             obj.cost.computeFunctionAndGradient();
%             obj.constraint.computeFunctionAndGradient();
%             Dg = obj.constraint.gradient;
%             S  = (Dg'*Dg)^-1;
% 
%             obj.designVariable.update(x);
%             obj.cost.computeFunctionAndGradient();
%             obj.constraint.computeFunctionAndGradient();
%             gx = obj.constraint.value;
%             lGx = max(-obj.p*obj.dualUpdater.lJPl,min(obj.p*obj.dualUpdater.lJPl,0.5*obj.aG*S*gx));
%             mF = mF+gx*(lGx-obj.dualUpdater.lGPl);
%         end

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

%            if max(abs(obj.meritGradient)) < obj.tol && obj.checkConstraint()
%                obj.hasConverged = true;
%            end


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