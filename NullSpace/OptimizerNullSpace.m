classdef OptimizerNullSpace < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
    end

    properties (Access = private)
        lineSearchTrials
        tol = 1e-8
        hasConverged
        acceptableStep
        hasFinished
        mOld
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
        gJFlowRatio
        predictedTau
        firstEstimation
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
            obj.hasFinished  = false;
            obj.printOptimizerVariable();
            obj.updateMonitoring();
            obj.computeNullSpaceFlow();
            obj.computeRangeSpaceFlow();
            obj.firstEstimation = false;
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
            obj.maxIter        = cParams.maxIter;
            obj.eta            = 0;
            obj.lG             = 0;
            obj.lJ             = 0;
            obj.etaNorm        = cParams.etaNorm;

            % Important parameters %%%%
            obj.etaMin         = 1e-6; % Just in case the null space flow is too small. May be analogous with Florian's n0 parameter in the future
            obj.etaMax         = cParams.etaMax; % For density can be Inf; for level-set this is adjusted with the density case, but we can increase it a bit
            obj.etaMaxMin      = cParams.etaMaxMin; % 'only TUNING' for level-set, in the end etaMax decreases until etaMaxMin; for density this 'only TUNING' is equivalent to tauMax
            obj.gJFlowRatio    = cParams.gJFlowRatio; % Robust parameter
            %%%%%%

            obj.hasConverged   = false;
            obj.nIter          = 0;
            obj.meritOld       = 1e6;
            obj.firstEstimation = true;
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
            titles  = [titles;{'Line Search';'Line Search trials';'Eta';'EtaMax';'lG';'lJ';'1/\eta (1-gk1/gk)';'Merit'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'log'},chConstr,{'bar','bar','plot','log','plot','plot','plot','plot'}];
            switch class(obj.designVariable)
                case 'LevelSet'
                    titles = [titles;{'Theta';'Alpha';'Beta';'Mean constr'}];
                    chartTypes = [chartTypes,{'plot','plot','plot','plot'}];
            end
            s.shallDisplay = cParams.monitoring;
            s.maxNColumns  = 6;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            data = obj.cost.value;
            data = [data;obj.cost.getFields(':')];
            data = [data;obj.constraint.value];
            data = [data;obj.designVariable.computeL2normIncrement()];
            data = [data;obj.dualVariable.fun.fValues];
            if obj.nIter == 0
                data = [data;0;0;0;obj.etaMax;0;0;0;NaN];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials;obj.eta;obj.etaMax;norm(obj.lG);norm(obj.lJ);norm(obj.predictedTau);obj.meritNew];
            end
            switch class(obj.designVariable)
                case 'LevelSet'
                    if obj.nIter == 0
                        data = [data;0;0;0;obj.constraint.gMean];
                    else
                        data = [data;obj.primalUpdater.Theta;obj.primalUpdater.Alpha;obj.primalUpdater.Beta;obj.constraint.gMean];
                    end
            end
            obj.monitoring.update(obj.nIter,data);
        end

        function updateEtaParameter(obj)
            vgJ     = obj.gJFlowRatio;
            l2DxJ   = norm(obj.DxJ);
            l2Dxg   = norm(obj.Dxg);
            obj.eta = max(min(vgJ*l2DxJ/l2Dxg,obj.etaMax),obj.etaMin);
            obj.updateMonitoringMultipliers();
        end

        function updateMonitoringMultipliers(obj) % Used just to monitor
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
            g0 = obj.constraint.value;
            x0 = obj.designVariable.fun.fValues;
            obj.updateEtaParameter();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.dualUpdater.update(obj.eta,obj.primalUpdater);
            obj.mOld = obj.computeMeritFunction();
            obj.computeNullSpaceFlow();
            obj.computeRangeSpaceFlow();
            obj.computeMeritGradient();
            obj.calculateInitialStep();
            while ~obj.acceptableStep
                obj.updatePrimal();
                obj.checkStep(x0,g0);
            end
        end

        function calculateInitialStep(obj)
            x   = obj.designVariable;
            DmF = obj.meritGradient;
            if obj.nIter == 0
                factor = 1000;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                factor = 1.02;
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

        function checkStep(obj,x0,g0)
            mNew = obj.computeMeritFunction();
            x    = obj.designVariable.fun.fValues;
            g    = obj.constraint.value;
            etaN = obj.obtainTrustRegion();
            if mNew <= obj.mOld+1e-3  &&  norm(x-x0)/norm(x0) < etaN
                obj.predictedTau   = (1-g/g0)/obj.eta;
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
                obj.dualUpdater.updateOld();
                obj.updateEtaMax(g,g0);
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.predictedTau   = (1-g/g0)/obj.eta;
                obj.acceptableStep = true;
                obj.meritNew       = obj.mOld;
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function updateEtaMax(obj,g,g0)
            switch class(obj.primalUpdater)
                case 'SLERP'
%                     active = not(obj.checkIndividualConstraint());
%                     g = g.*active;
%                     if min(g.*g0)<-1e-10
%                         obj.etaMax  = obj.etaMax/2;
%                         obj.etaNorm = max(obj.etaNorm/1.02,0.001);
%                     elseif min(g.*g0)>1e-10
%                         obj.etaMax = obj.etaMax*1.2;
%                     end
%                     dPsi = obj.designVariable.computeIncrement();
%                     s.fValues = obj.meritGradient;
%                     s.mesh = obj.designVariable.fun.mesh;
%                     s.order = 'P1';
%                     TD = LagrangianFunction(s);
%                     TDn = TD.normalize('L2');
%                     s.operation = @(xV) dPsi.evaluate(xV)./TDn.evaluate(xV);
%                     tFun = DomainFunction(s);
%                     tNorm = Norm.computeL2(obj.designVariable.fun.mesh,tFun);
%                     obj.etaMax = sqrt(tNorm);

                    % if obj.constraint.gMean < 1e-3
                    %     obj.etaMax = max(obj.etaMax/1.2,obj.etaMaxMin);
                    % end

                    [actg,~] = obj.computeActiveConstraintsGradient();
                    isAlmostFeasible  = norm(actg) < 0.01;
                    isAlmostOptimal   = abs(obj.meritNew - obj.meritOld) < 0.001;
                    if isAlmostFeasible && isAlmostOptimal
                        obj.etaMax = max(obj.etaMax/1.05,obj.etaMaxMin);
                    end

                case 'HAMILTON-JACOBI'
                    obj.etaMax = Inf; % Not verified
                otherwise
                    t          = obj.primalUpdater.tau;
                    obj.etaMax = 1/t;
            end
        end

        function areAcceptable = checkIndividualConstraint(obj)
            for i = 1:length(obj.constraint.value)
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        areAcceptable(i,1) = abs(obj.constraint.value(i)) < obj.tol;
                    case {'INEQUALITY'}
                        areAcceptable(i,1) = obj.constraint.value(i) < -1e-2;
                end
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
            if abs(obj.meritNew - obj.meritOld) < obj.tol && obj.checkConstraint()
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