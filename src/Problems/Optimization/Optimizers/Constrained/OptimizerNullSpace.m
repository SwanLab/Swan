classdef OptimizerNullSpace < handle

    properties (Access = private)
        tolCost   = 1e-8
        tolConstr = 1e-6
    end

    properties (Access = private)
        cost
        constraint
        constraintCase
        designVariable
        dualVariable
        primalUpdater
        dualUpdater
        maxIter
        nIter
        monitoring
        lineSearchTrials
        hasConverged
        acceptableStep
        hasFinished
        mOldPrimal
        meritNew
        meritOld
        meritGradient
        DxJ
        Dxg
        eta
        etaMin
        etaMax
        etaMaxMin
        lG
        lJ
        etaNorm
        etaNormMin
        gJFlowRatio
        firstEstimation
    end

    methods (Access = public) 
        function obj = OptimizerNullSpace(cParams)
            obj.init(cParams);
            obj.createMonitoring(cParams);
            obj.prepareFirstIter();
        end

        function solveProblem(obj)
            obj.hasConverged = false;
            obj.hasFinished  = false;
            obj.plotVariable();
            obj.updateMonitoring();
            obj.computeNullSpaceFlow();
            obj.computeRangeSpaceFlow();
            obj.firstEstimation = false;
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.plotVariable();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.designVariable.updateOld();
            end
        end
    end

    methods(Access = private)
        function init(obj,cParams)
            obj.cost            = cParams.cost;
            obj.constraint      = cParams.constraint;
            obj.constraintCase  = cParams.constraintCase;
            obj.designVariable  = cParams.designVariable;
            obj.maxIter         = cParams.maxIter;
            obj.lG              = 0;
            obj.lJ              = 0;
            obj.gJFlowRatio     = cParams.gJFlowRatio;
            obj.hasConverged    = false;
            obj.nIter           = 0;
            obj.meritOld        = 1e6;
            obj.firstEstimation = true;
            obj.etaNorm         = cParams.etaNorm;
            obj.eta             = 0;
            obj.etaMin          = 1e-6;
            obj.primalUpdater   = cParams.primalUpdater;
            obj.createDualVariable();
            cParams.dualVariable = obj.dualVariable; % ESTO LO QUITARÉ CUANDO DUALVARIABLE DESAPAREZCA DEL DUAPUPDATERNULL..
            obj.dualUpdater     = DualUpdaterNullSpace(cParams);
            obj.initOtherParameters(cParams);
        end

        function createDualVariable(obj)
            s.nConstraints   = length(obj.constraintCase);
            obj.dualVariable = DualVariable(s);
        end

        function initOtherParameters(obj,cParams)
            switch class(obj.designVariable)
                case 'LevelSet'
                    obj.etaMax     = cParams.etaMax;
                    obj.etaMaxMin  = cParams.etaMaxMin;
                    obj.etaNormMin = cParams.etaNormMin;
                otherwise
                    obj.etaMax = inf;
            end
        end

        function createMonitoring(obj,cParams)
            s.shallDisplay   = cParams.monitoring;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.primalUpdater  = obj.primalUpdater;
            obj.monitoring   = MonitoringNullSpace(s);
        end

        function updateMonitoring(obj)
            s.etaMax           = obj.etaMax;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.eta              = obj.eta;
            s.lG               = obj.lG;
            s.lJ               = obj.lJ;
            s.meritNew         = obj.meritNew;
            obj.monitoring.update(obj.nIter,s);
            obj.monitoring.refresh();
        end

        function plotVariable(obj)
            if ismethod(obj.designVariable,'plot')
                obj.designVariable.plot();
            end
        end

        function updateEtaParameter(obj)
            vgJ     = obj.gJFlowRatio;
            l2DxJ   = norm(obj.DxJ);
            l2Dxg   = norm(obj.Dxg);
            obj.eta = max(min(vgJ*l2DxJ/l2Dxg,obj.etaMax),obj.etaMin);
            obj.updateMonitoringMultipliers();
        end

        function updateMonitoringMultipliers(obj)
            g      = obj.constraint.value;
            Dg     = obj.constraint.gradient;
            DJ     = obj.cost.gradient;
            obj.lG = obj.eta*((Dg'*Dg)\g);
            obj.lJ = -1*((Dg'*Dg)\Dg')*DJ;
        end

        function computeNullSpaceFlow(obj)
            DJ     = obj.cost.gradient;
            [~,Dg] = obj.computeActiveConstraintsGradient();
            if isempty(Dg)
                obj.DxJ = DJ;
            else
                obj.DxJ = DJ-(Dg*(((Dg'*Dg)\Dg')*DJ));
            end
        end

        function computeRangeSpaceFlow(obj)
            [g,Dg] = obj.computeActiveConstraintsGradient();
            if isempty(Dg)
                obj.Dxg = zeros(size(obj.DxJ));
            else
                obj.Dxg = Dg*((Dg'*Dg)\g);
            end
        end

        function [actg,actDg] = computeActiveConstraintsGradient(obj)
            l   = obj.dualVariable.fun.fValues;
            gCases = obj.constraintCase;
            active = false(length(gCases),1);
            for i = 1:length(gCases)
                switch gCases{i}
                    case 'EQUALITY'
                        active(i) = 1;
                    case 'INEQUALITY'
                        if l(i)>1e-6 || obj.firstEstimation
                            active(i) = 1;
                        end
                end
            end
            actg   = obj.constraint.value(active);
            actDg  = obj.constraint.gradient(:,active);
        end

        function prepareFirstIter(obj)
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            obj.designVariable.updateOld();
        end

        function update(obj)
            x0 = obj.designVariable.fun.fValues;
            obj.updateEtaParameter();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.dualUpdater.update(obj.eta,obj.primalUpdater);
            obj.mOldPrimal = obj.computeMeritFunction();
            obj.computeNullSpaceFlow();
            obj.computeRangeSpaceFlow();
            obj.computeMeritGradient();
            obj.calculateInitialStep();
            while ~obj.acceptableStep
                obj.updatePrimal();
                obj.checkStep(x0);
            end
        end

        function calculateInitialStep(obj)
            x   = obj.designVariable;
            DmF = obj.meritGradient;
            if obj.nIter == 0
                factor = 50;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                factor = 1.05;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function updatePrimal(obj)
            x = obj.designVariable;
            g = obj.meritGradient;
            x = obj.primalUpdater.update(g,x);
            obj.designVariable = x;
        end

        function computeMeritGradient(obj)
            DJ  = obj.cost.gradient;
            Dg  = obj.constraint.gradient;
            l   = obj.dualVariable.fun.fValues;
            DmF = DJ+Dg*l;
            obj.meritGradient = DmF;
        end

        function checkStep(obj,x0)
            mNew = obj.computeMeritFunction();
            x    = obj.designVariable.fun.fValues;
            etaN = obj.obtainTrustRegion();
            if mNew <= obj.mOldPrimal+1e-3  &&  norm(x-x0)/(norm(x0)+1) < etaN
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
                obj.dualUpdater.updateOld();
                obj.updateEtaMax();
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew       = obj.mOldPrimal;
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function updateEtaMax(obj)
            switch class(obj.primalUpdater)
                case 'SLERP'
                    [actg,~] = obj.computeActiveConstraintsGradient();
                    isAlmostFeasible  = norm(actg) < 0.01;
                    isAlmostOptimal   = obj.primalUpdater.Theta < 0.15;
                    if isAlmostFeasible && isAlmostOptimal
                        obj.etaMax  = max(obj.etaMax/1.05,obj.etaMaxMin);
                        obj.etaNorm = max(obj.etaNorm/1.1,obj.etaNormMin);
                    end
                case 'HAMILTON-JACOBI'
                    obj.etaMax = Inf; % Not verified
                otherwise
                    t          = obj.primalUpdater.tau;
                    obj.etaMax = 1/t;
            end
        end

        function etaN = obtainTrustRegion(obj)
            switch class(obj.designVariable)
                case {'LevelSet','MultiLevelSet'}
                    if obj.nIter == 0
                        etaN = inf;
                    else
                        etaN = obj.etaNorm;
                    end
                otherwise
                    etaN = obj.etaNorm;
            end
        end

        function mF = computeMeritFunction(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            obj.constraint.computeFunctionAndGradient(x);
            l  = obj.dualVariable.fun.fValues;
            J  = obj.cost.value;
            h  = obj.constraint.value;
            mF = J+l'*h;
        end

        function obj = checkConvergence(obj)
            value = obj.constraint.value;
            cases = obj.constraintCase;
            if abs(obj.meritNew - obj.meritOld) < obj.tolCost && Optimizer.checkConstraint(value,cases,obj.tolConstr)
                obj.hasConverged = true;
                if obj.primalUpdater.isTooSmall()
                    obj.primalUpdater.tau = 1;
                end
            end
            obj.meritOld = obj.meritNew;
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