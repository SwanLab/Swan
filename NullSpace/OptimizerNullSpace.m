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
        ub
        lb
        eta
        lG
        lJ
        etaMax
        etaNorm
        gJFlowRatio
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
                obj.designVariable.updateOld();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.ub             = cParams.ub;
            obj.lb             = cParams.lb;
            obj.maxIter        = cParams.maxIter;
            obj.eta            = 0;
            obj.lG             = 0;
            obj.lJ             = 0;
            obj.etaNorm        = cParams.etaNorm;
            obj.gJFlowRatio    = cParams.gJFlowRatio;
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
            titles  = [titles;{'Line Search';'Line Search trials';'Eta';'lG';'lJ'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'log'},chConstr,{'bar','bar','plot','plot','plot'}];
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
                data = [data;0;0;0;0;0];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials;obj.eta;obj.lG;obj.lJ];
            end
            % merit?
            obj.monitoring.update(obj.nIter,data);
        end

        function updateEtaParameter(obj)
            obj.etaMax = 1;

            g = obj.constraint.value;
            Dg = obj.constraint.gradient;
            DJ = obj.cost.gradient;
            if norm(g)<=1e-3
                %obj.gJFlowRatio = 0;
            end
            if obj.nIter>0
                tau = obj.primalUpdater.tau;
            else
                tau = 1e-9;
            end
            vgJ     = obj.gJFlowRatio;
            DxJ     = obj.computeNullSpaceFlow();
            Dxg     = obj.computeRangeSpaceFlow();
            obj.eta = vgJ*min(DxJ/Dxg,obj.etaMax); % 2*DxJ/(h*tau))
            obj.lG  = obj.eta/(Dg'*Dg)*g;
            obj.lJ  = -1/(Dg'*Dg)*Dg'*DJ;
        end

        function DxJ = computeNullSpaceFlow(obj)
            DJ     = obj.cost.gradient;
            Dg     = obj.constraint.gradient;
            Prange = Dg*((Dg'*Dg)\Dg');
            DxJ    = norm((eye(size(Prange))-Prange)*DJ);
        end

        function Dxg = computeRangeSpaceFlow(obj)
            g   = obj.constraint.value;
            Dg  = obj.constraint.gradient;
            Dxg = norm(Dg*((Dg'*Dg)\g));
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
            obj.updateEtaParameter();
            obj.acceptableStep      = false;
            obj.lineSearchTrials    = 0;
            obj.dualUpdater.update(obj.eta,obj.ub,obj.lb);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.computeMeritGradient();
            obj.calculateInitialStep();

            while ~obj.acceptableStep
                x = obj.updatePrimal();
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x);
        end

        function calculateInitialStep(obj)
            x   = obj.designVariable;
            DmF = obj.meritGradient;
            if obj.nIter == 0
                factor = 1000;
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

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);
            etaN = obj.obtainTrustRegion();
            if mNew < obj.mOld && norm(x-x0)/norm(x0) < etaN
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

        function etaN = obtainTrustRegion(obj)
            switch class(obj.designVariable)
                case 'LevelSet'
                    if obj.nIter == 0
                        etaN = inf;
                    else
                        etaN = obj.etaNorm;
                    end
                otherwise
                    etaN = obj.etaNorm;
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