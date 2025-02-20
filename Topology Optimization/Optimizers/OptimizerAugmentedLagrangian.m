classdef OptimizerAugmentedLagrangian < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'Augmented Lagrangian';
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        costOld
        tol = 1e-8
        nX
        nConstr
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        hasFinished
        mOld
        meritNew
        penalty
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
        
        function obj = OptimizerAugmentedLagrangian(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
            obj.prepareFirstIter();
        end

        function obj = solveProblem(obj)
            x = obj.designVariable;
            obj.hasConverged = false;
            obj.cost.computeFunctionAndGradient(x);
            obj.constraint.computeFunctionAndGradient(x);
            obj.hasFinished = 0;
            obj.printOptimizerVariable();
            obj.updateMonitoring();
            while ~obj.hasFinished
                obj.update();
                obj.printOptimizerVariable();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.nConstr        = length(cParams.constraintCase);
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.nX             = obj.designVariable.fun.nDofs;
            obj.maxIter        = cParams.maxIter;
            obj.hasConverged   = false;
            obj.nIter          = 0;
            obj.penalty        = cParams.rho;
            obj.createMonitoring(cParams);
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            titlesConst   = obj.constraint.getTitleFields();
            nSFCost       = length(titlesF);
            nSFConstraint = length(titlesConst);
            titles        = [{'Cost'};titlesF;titlesConst;{'Norm L2 x';'Penalty'}];
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
            chartTypes = [{'plot'},chCost,chConstr,{'log'},{'plot'},chConstr,{'plot','bar','bar'}];
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
            data = [data;obj.penalty];
            data = [data;obj.dualVariable.fun.fValues];
            if obj.nIter == 0
                data = [data;0;0];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials];
            end
            obj.monitoring.update(obj.nIter,data);
        end

        function prepareFirstIter(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.fun.fValues = zeros(obj.nConstr,1);
        end

        function obj = update(obj)
            x0 = obj.designVariable.fun.fValues;
            obj.designVariable.update(x0);
            obj.saveOldValues(x0);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.computeMeritGradient();
            while ~obj.acceptableStep
                x = obj.updatePrimal();
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x.fun.fValues);
        end

        function displayIter(obj,x)
            m = obj.designVariable.mesh;
            bm = m.createBoundaryMesh();
            s.backgroundMesh = m;
            s.boundaryMesh   = bm;
            um = UnfittedMesh(s);
            um.compute(x);
            figure()
            um.plot();
        end

        function obj = calculateInitialStep(obj)
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            x       = obj.designVariable;
            l       = obj.dualVariable.fun.fValues;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            g       = obj.constraint.value;
            p       = obj.penalty;
            DmF     = DJ + Dg*(l + p*g);
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                factor = 1.05;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function x = updatePrimal(obj)
            x   = obj.designVariable;
            g   = obj.meritGradient;
            x   = obj.primalUpdater.update(g,x);
        end

        function computeMeritGradient(obj)
            Dh    = obj.constraint.gradient;
            DJ    = obj.cost.gradient;
            l     = obj.dualVariable.fun.fValues;
            p     = obj.penalty;
            gPlus = obj.defineConstraintValue();
            g     = (DJ + Dh*(l + p*gPlus));
            obj.meritGradient = g;
        end

        function mF = computeMeritFunction(obj,x)
            obj.designVariable.update(x)
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            J      = obj.cost.value;
            gPlus  = obj.defineConstraintValue();
            l      = obj.dualVariable.fun.fValues;
            rho    = obj.penalty;
            mF     = J + l'*gPlus + 0.5*rho*(gPlus'*gPlus);
        end

        function c = defineConstraintValue(obj)
            c   = obj.constraint.value;
            l   = obj.dualVariable.fun.fValues;
            rho = obj.penalty;
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case 'EQUALITY'
                        
                    case 'INEQUALITY'
                        c(i) = max(c(i),-l/rho);
                end
            end
        end

        function etaN = obtainTrustRegion(obj)
            switch class(obj.designVariable)
                case 'LevelSet'
                    if obj.nIter == 0
                        etaN = inf;
                    else
                        etaN = 0.02;
                    end
                otherwise
                    etaN = 0.02;
            end
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x.fun.fValues);
            etaN = obj.obtainTrustRegion();
            if mNew <= obj.mOld+1e-3  &&  norm(x.fun.fValues-x0)/(norm(x0)+1) < etaN
                obj.acceptableStep = true;
                obj.dualUpdater.updatePenalty(obj.penalty);
                obj.dualUpdater.update();
                obj.meritNew = mNew;
            elseif obj.primalUpdater.isTooSmall()
                %warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew = mNew; % Provisional value
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
        end

        function obj = updateOldValues(obj,x)
            obj.designVariable.update(x);
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
        end

        function obj = checkConvergence(obj)
           if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
               obj.hasConverged = true;
           else
               
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

        function saveVariablesForAnalysis(obj)
            i                           = obj.nIter + 1;
            obj.globalCost(i)           = obj.cost.value;
            obj.globalConstraint(:,i)   = obj.constraint.value;
            obj.globalCostGradient(i)   = norm(obj.cost.gradient);
            obj.globalMerit(i)          = obj.meritNew;
            obj.globalLineSearch(i)     = obj.primalUpdater.tau;
            obj.globalDual(:,i)         = obj.dualVariable.value;
            obj.globalDesignVar(:,i)    = obj.designVariable.value;
            if obj.hasConverged
                c = obj.globalCost;
                h = obj.globalConstraint;
                g = obj.globalCostGradient;
                m = obj.globalMerit;
                t = obj.globalLineSearch;
                d = obj.globalDual;
                v = obj.globalDesignVar;
                save('name.mat',"c","g","h","d","v");
            end
        end


    end

end