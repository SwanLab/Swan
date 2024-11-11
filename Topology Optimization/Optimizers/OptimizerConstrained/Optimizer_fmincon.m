classdef Optimizer_fmincon < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'fmincon';
    end

    properties (Access = private)
        problem
        upperBound
        lowerBound
        nX
        options
        algorithm
        hasConverged
        hasFinished
    end


    methods (Access = public)

        function obj = Optimizer_fmincon(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.createProblem();
            obj.createOptions();
        end

        function solveProblem(obj)
            f = obj.designVariable;
            obj.cost.computeFunctionAndGradient(f);
            obj.constraint.computeFunctionAndGradient(f);
            obj.designVariable.updateOld();
            obj.printOptimizerVariable();
            obj.updateMonitoring();
            x = obj.callfmincon();
            obj.designVariable.update(x);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.algorithm              = 'interior-point';
            cParams.optimizerNames.alg = obj.algorithm;
            obj.upperBound             = cParams.ub;
            obj.lowerBound             = cParams.lb;
            obj.nX                     = length(obj.designVariable.fun.fValues);
            obj.hasConverged           = false;
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
                chConstr{i}   = 'plot';
            end
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,chConstr,{'logy'}];

            s.shallDisplay = cParams.monitoring;
            s.maxNColumns  = 5;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            data = {};
            data{end+1} = obj.cost.value;
            data{end+1} = obj.cost.getFields(':');
            for i=1:length(obj.constraint.value)
                data{end+1} = obj.constraint.value(i);
            end
            data{end+1} = obj.designVariable.computeL2normIncrement();
            obj.monitoring.update(obj.nIter,data);
            obj.monitoring.refresh();
        end

        function x = callfmincon(obj)
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            x = fmincon(PROBLEM);
        end

        function createProblem(obj)
            prob.objective         = @(x) obj.objectiveAndGradient(x);
            prob.x0                = obj.designVariable.fun.fValues;
            prob.A                 = [];
            prob.b                 = [];
            prob.Aeq               = [];
            prob.beq               = [];
            prob.ub                = obj.upperBound*ones(obj.nX,1);
            prob.lb                = obj.lowerBound*ones(obj.nX,1);
            prob.nonlcon           = @(x) obj.constraintAndGradient(x);
            prob.solver            = "fmincon";
            obj.problem            = prob;
        end

        function [f,g] = objectiveAndGradient(obj,x)
            f = obj.objective(x);
            g = obj.gradient(x);
        end

        function f = objective(obj,x)
            obj.designVariable.update(x);
            f = obj.designVariable;
            obj.cost.computeFunctionAndGradient(f);
            f = obj.cost.value;
        end
        
        function g = gradient(obj,x)
            obj.designVariable.update(x);
            f = obj.designVariable;
            obj.cost.computeFunctionAndGradient(f);
            g = obj.cost.gradient;
        end

        function [c,ceq,DC,DCeq] = constraintAndGradient(obj,x)
            ceq    = [];
            DCeq   = [];
            c      = obj.constraintFunction(x);
            DC     = obj.constraint_gradient(x);
        end

        function f = constraintFunction(obj,x)
            obj.designVariable.update(x);
            f = obj.designVariable;
            obj.constraint.computeFunctionAndGradient(f)
            f = obj.constraint.value;
        end
        
        function g = constraint_gradient(obj,x)
            obj.designVariable.update(x);
            f = obj.designVariable;
            obj.constraint.computeFunctionAndGradient(f)
            g = obj.constraint.gradient;
        end

        function obj = createOptions(obj)
            obj.maxIter                    = 1e3;
            opts                           = optimoptions("fmincon");
            opts.Algorithm                 = obj.algorithm;
            opts.BarrierParamUpdate        = "monotone";
            opts.SpecifyObjectiveGradient  = true;
            opts.SpecifyConstraintGradient = true;
            opts.CheckGradients            = false;
            opts.ConstraintTolerance       = 1e-4;
            opts.EnableFeasibilityMode     = false;
            opts.HessianApproximation      = 'bfgs';
            opts.HessianFcn                = [];
            opts.HessianMultiplyFcn        = [];
            opts.HonorBounds               = true;
            opts.MaxFunctionEvaluations    = 100e3;
            opts.MaxIterations             = obj.maxIter;
            opts.StepTolerance	           = 1e-15;
            opts.OutputFcn                 = @(x,optimvalues,state)obj.myoutput(x,optimvalues,state);
            obj.options                    = opts;
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

    methods (Access = private)

        function stop = myoutput(obj,x,params,state)
            stop      = false;
            switch state
                case "init"

                case "iter"
                    obj.updateIterInfo();
                    params.algorithm   = obj.algorithm;
                    params.nIter       = obj.nIter;
                    params.hasFinished = obj.hasFinished; 
                    obj.printOptimizerVariable();
                    obj.updateMonitoring();
                    obj.designVariable.updateOld();
            end
        end

    end

end