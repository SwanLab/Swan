classdef OptimizerNullSpace < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
    end

    properties (Access = private)
        lineSearchTrials
        tol = 1e-5
        hasConverged
        acceptableStep
        hasFinished
        mOld
        meritNew
        meritGradient
        eta = 1
        etaNorm
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
            obj.prepareFirstIter();
        end

        function solveProblem(obj)
            obj.hasConverged = false;
            obj.hasFinished = false;
            obj.printOptimizerVariable();
            obj.updateMonitoring();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.printOptimizerVariable();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.checkParameters();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.maxIter        = cParams.maxIter;
            obj.etaNorm        = cParams.etaNorm;
            obj.hasConverged   = false;
            obj.nIter          = 0;
            obj.createMonitoring(cParams);
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            titlesConst   = obj.constraint.getTitleFields();
            nSFCost       = length(titlesF);
            nSFConstraint = length(titlesConst);
            titles        = [{'Cost'};titlesF;titlesConst;{'Norm L2 x'}];
            chConstr      = cell(1,nSFConstraint);
            for i = 1:nSFConstraint
                titles{end+1} = ['\lambda_{',titlesConst{i},'}'];
                chConstr{i}   = 'plot';
            end
            titles  = [titles;{'Line Search';'Line Search trials'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'log'},chConstr,{'bar','bar'}];
            s.shallDisplay = cParams.monitoring;
            s.maxNColumns  = 5;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            data = obj.cost.value;
            data = [data;obj.cost.getFields(':')];
            data = [data;obj.constraint.value];
            data = [data;obj.designVariable.computeL2normIncrement()];
            data = [data;obj.dualVariable.value];
            if obj.nIter == 0
                data = [data;0;0];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials];
            end
            % merit?
            obj.monitoring.update(obj.nIter,data);
        end

        function aJmax = nullSpaceParameterEstimation(obj,cParams)
            if isfield (cParams,'aJmax')
                aJmax = cParams.aJmax;
            else
                DJ = obj.cost.gradient;
                Dg = obj.constraint.gradient;
                aJmax = abs(-1/((Dg'*Dg)\Dg'*DJ));
            end
        end

        function checkParameters(obj)
            % recompute eta
        end

        function prepareFirstIter(obj)
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            obj.designVariable.updateOld();
            obj.dualVariable.value = zeros(size(obj.dualVariable.value));
        end

        function update(obj)
            x0 = obj.designVariable.fun.fValues;
            obj.calculateInitialStep();
            obj.acceptableStep      = false;
            obj.lineSearchTrials    = 0;
            obj.dualUpdater.update(obj.eta);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.computeMeritGradient();

            while ~obj.acceptableStep
                x = obj.updatePrimal();
                s.x  = x;
                s.x0 = x0;
                obj.checkStep(s);
            end
            obj.updateOldValues(x);
        end

        function calculateInitialStep(obj)
            x  = obj.designVariable.fun.fValues;
            DJ = obj.cost.gradient;
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(DJ,x,factor);
            else
                factor = 1.2;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function x = updatePrimal(obj)
            x       = obj.designVariable.fun.fValues;
            g       = obj.meritGradient;
            x       = obj.primalUpdater.update(g,x);
        end

        function computeMeritGradient(obj)
            DJ   = obj.cost.gradient;
            Dg   = obj.constraint.gradient;
            l   = obj.dualVariable.value;
            DmF = DJ+Dg*l;
            obj.meritGradient = DmF;
        end

        function checkStep(obj,s)
            x    = s.x;
            x0   = s.x0;
            mNew = obj.computeMeritFunction(x);
            if mNew < obj.mOld && norm(x-x0)/norm(x0) < obj.etaNorm
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
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

        function mF = computeMeritFunction(obj,xVal)
            x = obj.designVariable;
            x.update(xVal);
            obj.cost.computeFunctionAndGradient(x);
            obj.constraint.computeFunctionAndGradient(x);
            l  = obj.dualVariable.value;
            J  = obj.cost.value;
            h  = obj.constraint.value;
            mF = J+l'*h;
        end

        function obj = updateOldValues(obj,xV)
            x = obj.designVariable;
            x.update(xV);
            obj.cost.computeFunctionAndGradient(x);
            obj.constraint.computeFunctionAndGradient(x);
        end

        function obj = checkConvergence(obj)
            if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
                obj.hasConverged = true;
                if obj.primalUpdater.isTooSmall()
                    obj.primalUpdater.tau = 1;
                end
            end
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
            itHas = obj.nIter >= obj.maxIter;
        end

    end

end