classdef OptimizationProblem < handle

    properties (GetAccess = public, SetAccess = public)
        cost
        designVariable
        constraint
        dualVariable
        optimizer
    end

    methods (Access = public)

        function obj = OptimizationProblem(cParams)
            obj.init(cParams);
            obj.createOptimizer(cParams);
        end

        function solve(obj)

        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.designVariable = cParams.designVariable;
            obj.constraint     = cParams.constraint;
            obj.dualVariable   = cParams.dualVariable;
        end

        function createOptimizer(obj,cParams)
            % ...
            % also look what to do with Incremental scheme (I think it should be inside TopOpt_Problem)
        end

    end
end