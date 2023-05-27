classdef Volume_constraintWithBound < ShapeFunctional

properties (Access = public)
    volume
    vTarget
end

    methods (Access = public)

        function obj = Volume_constraintWithBound(cParams)
            cParams.designVariable = cParams.designVariable.density;
            obj.volume = Volume_constraint(cParams);
            obj.vTarget = cParams.targetParameters.Vfrac;
        end

        function computeFunctionAndGradient(obj)
            obj.volume.computeFunctionAndGradient();
            obj.value = obj.volume.value;
            obj.gradient = [obj.volume.gradient;0];
        end

        function v = getVariablesToPlot(obj)
            v{1} = (obj.volume.value+1)*obj.vTarget;
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'Volum';
        end

    end

end