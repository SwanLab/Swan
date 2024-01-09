classdef AcademicProblem < handle

    properties (Access = public)
        result
    end
    
    properties (Access = private)
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
            obj.solve();
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.createDesignVariable(cParams);
            obj.createCost(cParams);
            obj.createConstraint(cParams);
            obj.createOptimizer(cParams);
        end

        function createDesignVariable(obj,cParams)
            p.x0               = cParams.initialGuess;
            obj.designVariable = DesignVariableAcademic(p);
        end

        function createCost(obj,cParams)
            s.cF     = cParams.cost.cF;
            s.gF     = cParams.cost.gF;
            s.dV     = obj.designVariable;
            obj.cost = AcademicCost(s);
        end

        function createConstraint(obj,cParams)
            s.cF           = cParams.constraint.cF;
            s.gF           = cParams.constraint.gF;
            s.nSF          = cParams.constraint.nSF;
            s.dV           = obj.designVariable;
            obj.constraint = AcademicConstraint(s);
        end

        function createOptimizer(obj,cParams)
            s                            = cParams.settings;
            s.designVar                  = obj.designVariable;
            s.cost                       = obj.cost;
            s.constraint                 = obj.constraint;
            s.nConstraints               = obj.constraint.nSF;
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

    end

end