classdef LinearBoundFunction < ShapeFunctional

    methods (Access = public)

        function obj = LinearBoundFunction(cParams)
            obj.designVariable = cParams.designVariable;
        end

        function computeFunctionAndGradient(obj)
            obj.value         = obj.designVariable.value(end);
            obj.gradient      = zeros(length(obj.designVariable.value),1);
            obj.gradient(end) = 1;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Bound variable';
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value(end);
        end

    end

end