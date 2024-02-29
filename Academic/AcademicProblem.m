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
            obj.createAcademicCost(cParams);
            obj.createAcademicConstraint(cParams);
            obj.createOptimizer(cParams);
        end

        function createDesignVariable(obj,cParams)
            p.x0               = cParams.initialGuess;
            obj.designVariable = DesignVariableAcademic(p);
        end

        function createAcademicCost(obj,cParams)
            s.cF     = cParams.cost.cF;
            s.gF     = cParams.cost.gF;
            obj.cost = AcademicCost(s);
        end

        function createAcademicConstraint(obj,cParams)
            constCell = cell(length(cParams.constraint.cF),1);
            for i = 1:length(constCell)
                s.cF         = cParams.constraint.cF{i};
                s.gF         = cParams.constraint.gF{i};
                constCell{i} = AcademicConstraint(s);
            end
            obj.constraint = constCell;
        end

        function cost = createCost(obj)
            s.shapeFunctions{1} = obj.cost;
            s.weights           = 1;
            s.Msmooth           = 1;
            cost                = Cost(s);
        end

        function constraint = createConstraint(obj)
            s.shapeFunctions = obj.constraint;
            s.Msmooth        = 1;
            constraint       = Constraint(s);
        end

        function createOptimizer(obj,cParams)
            s                = cParams.settings;
            s.designVariable = obj.designVariable;
            s.cost           = obj.createCost();
            s.constraint     = obj.createConstraint();
            s.nConstraints   = length(cParams.constraint.cF);
            s.dualVariable   = DualVariable(s);
            s.monitoring     = true;
            s.primal         = 'PROJECTED GRADIENT';
            s.tolerance      = 1e-8;
            obj.optimizer    = Optimizer.create(s);
        end

        function solve(obj)
            obj.optimizer.solveProblem();
            obj.result = obj.designVariable.value;
        end

    end

end