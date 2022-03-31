classdef Optimizer_fmincon < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'fmincon';
    end

    properties
        problem
        nConstr
        info
        constraintTolerance
        optimalityTolerance
        data
        upperBound
        lowerBound
        functions
        nX
        options
        algorithm
    end

    methods (Access = public)

        function obj = Optimizer_fmincon(cParams)
            obj.algorithm  = 'interior-point';
            cParams.monitoringDockerSettings.optimizerNames.alg = obj.algorithm;
            cParams.optimizerNames.alg = obj.algorithm;
            obj.init(cParams);
            obj.upperBound = cParams.uncOptimizerSettings.ub;
            obj.lowerBound = cParams.uncOptimizerSettings.lb;
            obj.nConstr    = cParams.nConstr;
            obj.maxIter    = cParams.maxIter;
            obj.nIter      = -1;
            obj.nX         = length(obj.designVariable.value);
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

        function x = callfmincon(obj)
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            x = fmincon(PROBLEM);
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
            opts.Display                   = "none";
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

    end

    methods (Access = private)

        function stop = myoutput(obj,x,params,state)
            stop      = false;
            switch state
                case "init"
                    
                case "iter"
                    obj.data  = params;
                    obj.nIter = obj.nIter+1;
                    obj.designVariable.update(x);
                    foOpt       = params.firstorderopt;
                    normXsquare = obj.designVariable.computeL2normIncrement();
                    obj.designVariable.updateOld();
                    incX = sqrt(normXsquare);

                    switch obj.algorithm
                        case 'sqp'
                            stepL = params.stepsize;
                        case 'interior-point'
                            stepL = params.trustregionradius;
                        otherwise
                    end

                    obj.updateStatus();
                    obj.printOptimizerVariable();
                    obj.dualVariable.value = zeros(obj.constraint.nSF,1);            
                    obj.convergenceVars.reset();
                    obj.convergenceVars.append(incX);
                    obj.convergenceVars.append(foOpt);
                    obj.convergenceVars.append(stepL);
                    obj.refreshMonitoring();
                    obj.printHistory();
                otherwise
                    
            end
        end

    end



end