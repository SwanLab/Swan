classdef Optimizer_fmincon < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'fmincon';
    end

    properties (Access = private)
        problem
        iterDisplay
        upperBound
        lowerBound
        nX
        options
        algorithm
        incrementalScheme
        hasConverged
        hasFinished
        
        globalCost
        globalConstraint
        globalCostGradient
        globalLineSearch
        globalDual
        globalDesignVar
    end


    methods (Access = public)

        function obj = Optimizer_fmincon(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createProblem();
            obj.createOptions();
        end

         function solveProblem(obj)
            obj.cost.computeFunctionAndGradient();
            obj.designVariable.updateOld();                
            x = obj.callfmincon();
            obj.designVariable.update(x);
         end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.algorithm              = 'interior-point';
            cParams.optimizerNames.alg = obj.algorithm;
            obj.upperBound             = cParams.uncOptimizerSettings.ub;
            obj.lowerBound             = cParams.uncOptimizerSettings.lb;
            obj.iterDisplay            = cParams.outputFunction.iterDisplay;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.nX                     = length(obj.designVariable.value);
            obj.hasConverged           = false;
            cParams.monitoringDockerSettings.optimizerNames.alg = obj.algorithm;
        end

        function x = callfmincon(obj)
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            x = fmincon(PROBLEM);
            v = obj.globalDesignVar;
            c = obj.globalCost;
            h = obj.globalConstraint;
            save('fminconIPOPTacademic4','v','c','h');
        end

        function createProblem(obj)
            prob.objective         = @(x) obj.objectiveAndGradient(x);
            prob.x0                = obj.designVariable.value;
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
            obj.cost.computeFunctionAndGradient();
            f = obj.cost.value;
        end
        
        function g = gradient(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient()
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
            obj.constraint.computeFunctionAndGradient()
            f = obj.constraint.value;
        end
        
        function g = constraint_gradient(obj,x)
            obj.designVariable.update(x);
            obj.constraint.computeFunctionAndGradient()
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
            opts.Display                   = obj.iterDisplay;
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
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end

    end

    methods (Access = private)

        function stop = myoutput(obj,x,params,state)
            stop      = false;
            switch state
                case "init"
                    obj.globalDesignVar(:,obj.nIter + 1) = x;
                case "iter"
                    obj.updateIterInfo();
                    obj.designVariable.update(x);
                    obj.globalDesignVar(:,obj.nIter + 1) = x;
                    obj.globalCost(obj.nIter+1)       = obj.cost.value;
                    obj.globalConstraint(:,obj.nIter+1) = obj.constraint.value;
                    params.algorithm   = obj.algorithm;
                    params.nIter       = obj.nIter;
                    params.hasFinished = obj.hasFinished; 
                    obj.outputFunction.monitoring.compute(params);
            end
        end

    end

end