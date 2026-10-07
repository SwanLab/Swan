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
        gif
        gifName
        printing
        printName
        compliance

    end

    properties (Access = public)

        costHistory       = [];
        complianceHistory = [];
        volumeConstraintHistory = [];
        tauHistory              = [];
        etaHistory              = [];
        etaMaxHistory           = [];
        lineSearchTrialsHistory = [];

        % =========================================================
        % Two-variable optimization diagnostics
        % =========================================================
        normDJbHistory       = [];
        normDJrhoHistory     = [];

        normDxJbHistory      = [];
        normDxJrhoHistory    = [];

        normMeritBHistory    = [];
        normMeritRhoHistory  = [];

        deltaBHistory        = [];
        deltaRhoHistory      = [];

        fracBLowerHistory    = [];
        fracBUpperHistory    = [];

        fracRhoLowerHistory  = [];
        fracRhoUpperHistory  = [];
        fracRhoBoundHistory  = [];

        meanBHistory         = [];
        meanRhoHistory       = [];

        minBHistory          = [];
        maxBHistory          = [];

        minRhoHistory        = [];
        maxRhoHistory        = [];
        meritConditionHistory = [];
        trustConditionHistory = [];

        meritDifferenceHistory = [];
        trustRatioHistory = [];

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
                obj.printResults();
                obj.updateIterInfo();
                obj.plotVariable();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.designVariable.updateOld();
            end
        end
        function data = getFinalKKTDiagnostics(obj)

            n = length(obj.designVariable.funB.fValues);

            x = obj.designVariable.fun.fValues;
            g = obj.meritGradient;

            lb = [-ones(n,1); 1e-6*ones(n,1)];
            ub = [ ones(n,1); 0.998*ones(n,1)];

            % Unit projected-gradient test
            alpha = 1;

            xProj = min(ub,max(x-alpha*g,lb));

            r = x-xProj;

            data.rB   = r(1:n);
            data.rRho = r(n+1:2*n);

            data.normRB   = norm(data.rB);
            data.normRRho = norm(data.rRho);
            data.normR    = norm(r);

            data.maxRB   = max(abs(data.rB));
            data.maxRRho = max(abs(data.rRho));

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
            obj.gif             = cParams.gif;
            obj.gifName         = cParams.gifName;
            obj.printing        = cParams.printing;
            obj.printName       = cParams.printName;
            obj.primalUpdater   = cParams.primalUpdater;
            obj.dualUpdater     = DualUpdaterNullSpace(cParams);
            obj.createDualVariable();
            obj.initOtherParameters(cParams);
            if isfield(cParams,'compliance')
                obj.compliance = cParams.compliance;
            else
                obj.compliance = [];
            end


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

        % function updateMonitoring(obj)
        %     s.etaMax           = obj.etaMax;
        %     s.lineSearchTrials = obj.lineSearchTrials;
        %     s.eta              = obj.eta;
        %     s.lG               = obj.lG;
        %     s.lJ               = obj.lJ;
        %     s.meritNew         = obj.meritNew;
        %     obj.monitoring.update(obj.nIter,s);
        %     obj.monitoring.refresh();
        % end
        function updateMonitoring(obj)

            s.etaMax           = obj.etaMax;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.eta              = obj.eta;
            s.lG               = obj.lG;
            s.lJ               = obj.lJ;
            s.meritNew         = obj.meritNew;

            obj.monitoring.update(obj.nIter,s);
            obj.monitoring.refresh();

            costValue = obj.cost.value;
            costValue = costValue(1);

            obj.costHistory(end+1,1) = costValue;
            obj.complianceHistory(end+1,1) = costValue;

            constraintValue = obj.constraint.value;
            constraintValue = constraintValue(1);

            obj.volumeConstraintHistory(end+1,1) = constraintValue;

            tauNow = obj.primalUpdater.tau;

            if isempty(tauNow)
                obj.tauHistory(end+1,1) = NaN;
            else
                obj.tauHistory(end+1,1) = tauNow;
            end

            if isempty(obj.eta)
                obj.etaHistory(end+1,1) = NaN;
            else
                obj.etaHistory(end+1,1) = obj.eta;
            end

            if isempty(obj.etaMax)
                obj.etaMaxHistory(end+1,1) = NaN;
            else
                obj.etaMaxHistory(end+1,1) = obj.etaMax;
            end

            if isempty(obj.lineSearchTrials)
                obj.lineSearchTrialsHistory(end+1,1) = NaN;
            else
                obj.lineSearchTrialsHistory(end+1,1) = ...
                    obj.lineSearchTrials;
            end

            hasTwoVariables = ...
                isprop(obj.designVariable,'funB') && ...
                isprop(obj.designVariable,'funRho');

            if hasTwoVariables

                n = length(obj.designVariable.funB.fValues);

                DJ = obj.cost.gradient;

                if ~isempty(DJ) && length(DJ) >= 2*n

                    DJb   = DJ(1:n);
                    DJrho = DJ(n+1:2*n);

                    obj.normDJbHistory(end+1,1)   = norm(DJb);
                    obj.normDJrhoHistory(end+1,1) = norm(DJrho);

                else

                    obj.normDJbHistory(end+1,1)   = NaN;
                    obj.normDJrhoHistory(end+1,1) = NaN;

                end

                if ~isempty(obj.DxJ) && length(obj.DxJ) >= 2*n

                    DxJb   = obj.DxJ(1:n);
                    DxJrho = obj.DxJ(n+1:2*n);

                    obj.normDxJbHistory(end+1,1)   = norm(DxJb);
                    obj.normDxJrhoHistory(end+1,1) = norm(DxJrho);

                else

                    obj.normDxJbHistory(end+1,1)   = NaN;
                    obj.normDxJrhoHistory(end+1,1) = NaN;

                end

                Dm = obj.meritGradient;

                if ~isempty(Dm) && length(Dm) >= 2*n

                    DmB   = Dm(1:n);
                    DmRho = Dm(n+1:2*n);

                    obj.normMeritBHistory(end+1,1)   = norm(DmB);
                    obj.normMeritRhoHistory(end+1,1) = norm(DmRho);

                else

                    obj.normMeritBHistory(end+1,1)   = NaN;
                    obj.normMeritRhoHistory(end+1,1) = NaN;

                end

            else

                obj.normDJbHistory(end+1,1)       = NaN;
                obj.normDJrhoHistory(end+1,1)     = NaN;

                obj.normDxJbHistory(end+1,1)      = NaN;
                obj.normDxJrhoHistory(end+1,1)    = NaN;

                obj.normMeritBHistory(end+1,1)    = NaN;
                obj.normMeritRhoHistory(end+1,1)  = NaN;

            end

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
            obj.updateDualVariable();
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

        function printResults(obj)
            if obj.nIter/10==round(obj.nIter/10)
                if obj.gif
                    obtainGIF(obj.gifName,obj.designVariable,obj.nIter);
                end
                if obj.printing
                    obj.designVariable.fun.print([obj.printName,'Iter',num2str(obj.nIter/10)]);
                end
            end
        end

        function updateDualVariable(obj)
            if obj.nIter == 0
                lUB = 0;
                lLB = 0;
            else
                t   = obj.primalUpdater.boxConstraints.refTau;
                lUB = obj.primalUpdater.boxConstraints.lUB/t;
                lLB = obj.primalUpdater.boxConstraints.lLB/t;
            end
            l   = obj.dualUpdater.update(obj.eta,lUB,lLB);
            obj.dualVariable.update(l);
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

            % ========================================================
            % LINE SEARCH DIAGNOSTICS
            % ========================================================

            meritDiff = mNew - obj.mOldPrimal;

            trustRatio = ...
                norm(x-x0)/(norm(x0)+1);

            meritOK = meritDiff <= 1e-1;

            trustOK = trustRatio < etaN;


            % ========================================================
            % ACCEPT STEP
            % ========================================================

            if meritOK && trustOK

                obj.acceptableStep = true;

                obj.meritNew = mNew;

                obj.updateEtaMax();


                % ========================================================
                % STEP BECAME TOO SMALL
                % ========================================================

            elseif obj.primalUpdater.isTooSmall()

                warning( ...
                    'Convergence could not be achieved (step length too small)')

                fprintf('\n=================================================\n');
                fprintf('STEP TOO SMALL\n');
                fprintf('iter       = %d\n',obj.nIter);
                fprintf('trial      = %d\n',obj.lineSearchTrials);
                fprintf('tau        = %.6e\n',obj.primalUpdater.tau);
                fprintf('meritDiff  = %.6e\n',meritDiff);
                fprintf('meritOK    = %d\n',meritOK);
                fprintf('trustRatio = %.6e\n',trustRatio);
                fprintf('etaN       = %.6e\n',etaN);
                fprintf('trustOK    = %d\n',trustOK);
                fprintf('=================================================\n');

                obj.acceptableStep = true;

                obj.meritNew = obj.mOldPrimal;

                obj.designVariable.update(x0);


                % ========================================================
                % REJECT STEP -> REDUCE TAU
                % ========================================================

            else

                % Print only when many reductions are already happening
                if obj.lineSearchTrials >= 20

                    fprintf(['\nLS diagnostic: ' ...
                        'iter=%d trial=%d ' ...
                        'tau=%.3e ' ...
                        'meritDiff=%.3e meritOK=%d ' ...
                        'trustRatio=%.3e etaN=%.3e trustOK=%d\n'], ...
                        obj.nIter, ...
                        obj.lineSearchTrials, ...
                        obj.primalUpdater.tau, ...
                        meritDiff, ...
                        meritOK, ...
                        trustRatio, ...
                        etaN, ...
                        trustOK);

                end
                if obj.lineSearchTrials == 25

                    % Current candidate
                    xTrial = obj.designVariable.fun.fValues;

                    % Merit at trial point
                    meritTrial = mNew;

                    % --------------------------------------------------------
                    % Re-evaluate merit EXACTLY at x0
                    % --------------------------------------------------------
                    obj.designVariable.update(x0);

                    meritAtX0Again = obj.computeMeritFunction();
                    Jagain = obj.cost.value;
                    gagain = obj.constraint.value;
                    lambda = obj.dualVariable.fun.fValues;

                    fprintf('\nComponents at x0 re-evaluation:\n');
                    fprintf('J      = %.15e\n',Jagain(1));
                    fprintf('g      = %.15e\n',gagain(1));
                    fprintf('lambda = %.15e\n',lambda(1));
                    fprintf('lambda*g = %.15e\n',lambda(1)*gagain(1));

                    fprintf('\n=================================================\n');
                    fprintf('ZERO-STEP CONSISTENCY TEST\n');
                    fprintf('iter              = %d\n',obj.nIter);
                    fprintf('tau               = %.6e\n',obj.primalUpdater.tau);
                    fprintf('||xTrial-x0||     = %.6e\n',norm(xTrial-x0));
                    fprintf('\n');
                    fprintf('mOldPrimal        = %.15e\n',obj.mOldPrimal);
                    fprintf('meritTrial        = %.15e\n',meritTrial);
                    fprintf('meritAtX0Again    = %.15e\n',meritAtX0Again);
                    fprintf('\n');
                    fprintf('trial-old         = %.15e\n', ...
                        meritTrial-obj.mOldPrimal);
                    fprintf('x0again-old       = %.15e\n', ...
                        meritAtX0Again-obj.mOldPrimal);
                    fprintf('=================================================\n');

                    % Restore trial point so checkStep logic is not altered
                    obj.designVariable.update(xTrial);

                end

                obj.primalUpdater.decreaseStepLength();

                obj.designVariable.update(x0);

                obj.lineSearchTrials = ...
                    obj.lineSearchTrials + 1;

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