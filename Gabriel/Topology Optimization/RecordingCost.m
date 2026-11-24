classdef RecordingCost < handle

    properties (Dependent)
        value
        gradient
    end

    properties (Access = private)
        cost
        values
        fingerprints
    end

    methods

        function v = get.value(obj)
            v = obj.cost.value;
        end

        function g = get.gradient(obj)
            g = obj.cost.gradient;
        end

    end

    methods (Access = public)

        function obj = RecordingCost(cost)
            obj.cost         = cost;
            obj.values       = [];
            obj.fingerprints = zeros(0,4);
        end

        function computeFunctionAndGradient(obj,x)
            obj.cost.computeFunctionAndGradient(x);
            obj.record(x);
        end

        function nF = obtainNumberFields(obj)
            nF = obj.cost.obtainNumberFields();
        end

        function titles = getTitleFields(obj)
            titles = obj.cost.getTitleFields();
        end

        function j = getFields(obj,i)
            j = obj.cost.getFields(i);
        end

        function J = getEvaluationHistory(obj)
            J = obj.values;
        end

        function J = getIterationHistory(obj)
            f        = obj.fingerprints;
            repeated = all(f(2:end,:) == f(1:end-1,:),2);
            accepted = [repeated; true];
            J        = obj.values(accepted);
        end

    end

    methods (Access = private)

        function record(obj,x)
            xV = x.fun.fValues(:);
            obj.values(end+1,1)       = obj.cost.value;
            obj.fingerprints(end+1,:) = [sum(xV),xV'*xV,xV(1),xV(end)];
        end

    end

end