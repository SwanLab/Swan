classdef AcademicProblem < handle

    properties (Access = public)
        result
    end
    
    properties (Access = private)
        filename
        designVariable
        cost
        constraint
        optimizer
    end

    methods (Access = public)

        function obj = AcademicProblem(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.createOptimizer();
            obj.solve();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.filename = cParams.filename;
        end

        function createOptimizer(obj)
            pParams = obj.uploadProblemParameters();
            s       = pParams.settings;

            obj.prepareDesignVariable(pParams.initialGuess);
            obj.prepareCostConstraint(pParams);

            s.designVar                  = obj.designVariable;
            s.cost                       = obj.cost;
            s.constraint                 = obj.constraint;
            s.constraint.nSF             = pParams.nConstr;
            s.nConstraints               = pParams.nConstr;
            s.dualVariable               = DualVariable(s);
            s.outputFunction.type        = "Academic";
            s.outputFunction.iterDisplay = "iter";
            s.outputFunction.monitoring  = MonitoringManager(s);
            s.optimizerNames.primal     = 'PROJECTED GRADIENT';
            obj.optimizer = Optimizer.create(s);
        end

        function solve(obj)
            obj.optimizer.solveProblem();
            obj.result = obj.designVariable.value;
        end

        function pParams = uploadProblemParameters(obj)
            run(obj.filename);
            pParams.initialGuess     = x0;
            pParams.costHandle       = cost;
            pParams.constraintHandle = constraint;
            pParams.nConstr          = nConstr;
            pParams.settings         = s;
        end

        function prepareDesignVariable(obj,x0)
            p.x0               = x0;
            obj.designVariable = DesignVariableAcademic(p);
        end

        function prepareCostConstraint(obj,pParams)
            j              = pParams.costHandle;
            c              = pParams.constraintHandle;
            j.dV           = obj.designVariable;
            c.dV           = obj.designVariable;
            obj.cost       = AcademicCost(j);
            obj.constraint = AcademicConstraint(c);
        end

    end

end