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
        eta
        etaMax
        lG
        lJ
        etaNorm
        gJFlowRatio
        predictedTau
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
            obj.etaMax         = Inf;
            obj.etaNorm        = cParams.etaNorm;
            obj.gJFlowRatio    = cParams.gJFlowRatio;
            obj.hasConverged   = false;
            obj.nIter          = 0;
            obj.createMonitoring(cParams);
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            titlesConst   = obj.constraint.getTitleFields();
            titlesConst   = titlesConst{1,1};
            nSFCost       = length(titlesF);
            nSFConstraint = length(titlesConst);
            titles        = [{'Cost'};titlesF;titlesConst];
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
            chartTypes = [{'plot'},chCost,chConstr,chConstr,{'bar','bar','plot','plot','plot','plot','plot','plot'}];
            switch class(obj.designVariable)
                case 'LevelSet'
                    titles = [titles;{'Theta';'Alpha';'Beta'}];
                    chartTypes = [chartTypes,{'plot','plot','plot'}];
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
            data = [data;obj.dualVariable.fun.fValues];
            if obj.nIter == 0
                data = [data;0;0;0;obj.etaMax;0;0;0;NaN];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials;obj.eta;obj.etaMax;norm(obj.lG);norm(obj.lJ);norm(obj.predictedTau);obj.meritNew];
            end
            switch class(obj.designVariable)
                case 'LevelSet'
                    if obj.nIter == 0
                        data = [data;0;0;0];
                    else
                        data = [data;obj.primalUpdater.Theta;obj.primalUpdater.Alpha;obj.primalUpdater.Beta];
                    end
            end
            obj.monitoring.update(obj.nIter,data);
        end

        function updateEtaParameter(obj)
            vgJ     = obj.gJFlowRatio;
            DxJ     = obj.computeNullSpaceFlow();
            Dxg     = obj.computeRangeSpaceFlow();
            obj.eta = min(vgJ*DxJ/Dxg,obj.etaMax);
            obj.updateMonitoringMultipliers();
        end

        function updateMonitoringMultipliers(obj) % Used just to monitor
            g      = obj.constraint.value;
            Dg     = obj.constraint.gradient;
            DJ     = obj.cost.gradient;
            obj.lG = obj.eta*((Dg'*Dg)\g);
            obj.lJ = -1*((Dg'*Dg)\Dg')*DJ;
        end

        function DxJ = computeNullSpaceFlow(obj)
            DJ  = obj.cost.gradient;
            Dg  = obj.constraint.gradient;
            DxJ = norm(DJ-(Dg*(((Dg'*Dg)\Dg')*DJ)));
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
        end

        function update(obj)
            g0 = obj.constraint.value;
            x0 = obj.designVariable.fun.fValues;
            obj.updateEtaParameter();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.dualUpdater.update(obj.eta,obj.primalUpdater);
            obj.mOld = obj.computeMeritFunction();
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
                factor = 3;
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
            if mNew < obj.mOld && norm(x-x0)/norm(x0) < etaN
                obj.predictedTau   = (1-g/g0)/obj.eta;
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
                obj.dualUpdater.updateOld();
                obj.updateEtaMax();
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.predictedTau   = (1-g/g0)/obj.eta;
                obj.acceptableStep = true;
                obj.meritNew       = obj.mOld;
                x0 = reshape(x0,obj.designVariable.designVariable{1,1}.fun.mesh.nnodes,[]);
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
            else
                obj.primalUpdater.decreaseStepLength();
                x0 = reshape(x0,obj.designVariable.designVariable{1,1}.fun.mesh.nnodes,[]);
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function updateEtaMax(obj)
            switch class(obj.primalUpdater)
                case 'SLERP'
                    %                     phi = obj.designVariable.fun;
                    %                     quad = Quadrature.create(phi.mesh,'QUADRATIC');
                    %                     dV = phi.mesh.computeDvolume(quad);
                    %                     chiBall=obj.designVariable.computeBallCharacteristicFunction(quad);
                    %                     intBall = sum(dV.*chiBall,"all");
                    %                     gk = LagrangianFunction.create(phi.mesh,1,'P1');
                    %                     gk.fValues = obj.meritGradient;
                    %                     gkL2 = sqrt(Norm.computeL2(phi.mesh,gk));
                    %                     k  = obj.primalUpdater.tau;
                    %                     t  = obj.primalUpdater.Theta;
                    %                     obj.etaMax = intBall*gkL2*sin(t)/sin(k*t); % ak or not ak?
%                     theta      = obj.primalUpdater.Theta;
%                     k          = obj.primalUpdater.tau;
%                     b          = obj.primalUpdater.Beta;
%                     a          = obj.primalUpdater.Alpha;
                    obj.etaMax = Inf;
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