classdef OptimizerNullSpace < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        costOld
        tol = 1e-5
        nX
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        hasFinished
        mOld
        meritNew
        meritGradient
        aJmax
        aGmax
        aJ
        aG
        eta
        Vtar
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
            obj.prepareFirstIter();
            obj.aJmax = obj.nullSpaceParameterEstimation(cParams);
            obj.aGmax = obj.rangeSpaceParameterEstimation(cParams);
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
                %obj.updateMonitoring();
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
            obj.nX             = obj.designVariable.designVariable{1,1}.fun.nDofs;
            obj.maxIter        = cParams.maxIter;
            obj.hasConverged   = false;
            obj.nIter          = 0;
            obj.Vtar           = cParams.volumeTarget;
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
            titles  = [titles;{'Volume';'Line Search';'Line Search trials'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'log'},chConstr,{'plot','bar','bar'}];
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
            data = [data;obj.computeVolume(obj.constraint.value)]; % millorar
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

        function aGmax = rangeSpaceParameterEstimation(obj,cParams)
            if isfield (cParams,'aGmax')
                aGmax = cParams.aGmax;
            else
                Dg = obj.constraint.gradient;
                aGmax = 150/(inv(Dg'*Dg));
            end
        end

        function checkParameters(obj)
            if abs(obj.meritNew - obj.mOld) < 10*obj.tol
                g  = obj.constraint.value;
                DJ = obj.cost.gradient;
                if obj.aG <= 0.5*obj.aGmax
                    exponent = -sign(g)*sign(sum(DJ));
                    obj.aJmax = obj.aJmax*1^exponent;
                else
                    exponent  = 1-sign(obj.checkConstraint());
                    obj.aGmax = obj.aGmax*1^exponent;
                end
            end
        end

        function prepareFirstIter(obj)
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.value = zeros(size(obj.dualVariable.value));
        end

        function updateNullSpaceCoefficient(obj)
            targetVolume = obj.Vtar;
            obj.aJ       = obj.aJmax*(1-targetVolume);
        end

        function updateRangeSpaceCoefficient(obj)
            if obj.cost.value == 0
                obj.aG = obj.aGmax;
            else
                targetVolume = obj.Vtar;
                g            = obj.constraint.value;
                v            = obj.computeVolume(g); % class(obj.constraint.shapeFunctions{1,1})=='Volume_constraint'
                if v>targetVolume
                    r = 1-targetVolume;
                else
                    r = targetVolume;
                end
                obj.aG       = obj.aGmax*(1-abs(v-targetVolume)/r)^10;
            end
        end

        function updateMaximumVolumeRemoved(obj)
            if obj.nIter==0
                obj.eta = 0.05;
            else
                if obj.aG <= 0.5*obj.aGmax
                    obj.eta = 0.05;
                else
                    obj.eta = 0.01;
                end
            end
        end

        function update(obj)
            obj.updateNullSpaceCoefficient();
            obj.updateRangeSpaceCoefficient();
            obj.updateMaximumVolumeRemoved();
            x0 = obj.designVariable.fun.fValues;
            g0 = obj.constraint.value;
            obj.acceptableStep      = false;
            obj.lineSearchTrials    = 0;
            d.nullSpaceCoefficient  = obj.aJ;
            d.rangeSpaceCoefficient = obj.aG;
            obj.dualUpdater.update(d);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.computeMeritGradient();
            obj.calculateInitialStep();

            while ~obj.acceptableStep
                x = obj.updatePrimal();
                s.x  = x;
                s.x0 = x0;
                s.g0 = g0;
                obj.checkStep(s);
            end
            obj.updateOldValues(x);
        end

        function calculateInitialStep(obj)
            x   = obj.designVariable;
            DmF = obj.meritGradient;
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
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
            g0   = s.g0;
            mNew = obj.computeMeritFunction(x);
            g    = obj.constraint.value;
            v0   = obj.computeVolume(g0);
            v    = obj.computeVolume(g);
            if mNew < obj.mOld && norm(v-v0) < obj.eta
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

        function v = computeVolume(obj,g)
            targetVolume = obj.Vtar;
            v            = targetVolume*(1+g);
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