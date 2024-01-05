classdef DensityAndBound < DesignVariable

    properties (Access = public)
        density
        bound
    end

    methods (Access = public)

        function obj = DensityAndBound(cParams)
            obj.density = Density(cParams);
            obj.bound   = 1;
            obj.mesh    = obj.density.mesh;
            obj.type    = obj.density.type;
            obj.createValue();
        end

        function update(obj,value)
            obj.update@DesignVariable(value);
            obj.density.value = obj.value(1:end-1);
            obj.bound         = obj.value(end);
        end

        function updateOld(obj)
            obj.density.updateOld();
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.density.value;
        end

        function [fun, funNames] = getFunsToPlot(obj)
            [fun, funNames] = obj.density.getFunsToPlot();
        end

        function rho = computeVolumeFraction(obj)
            rho = obj.density.computeVolumeFraction();
        end

        function norm = computeL2normIncrement(obj)
           norm = obj.density.computeL2normIncrement();
        end

    end

    methods (Access = private)

        function createValue(obj)
            obj.value = [obj.density.fun.fValues;obj.bound]; % Habrá que crear una composite function FeAndBoundFun con fValues concatenados. Tmb la BoundFun por separado
        end
        
    end
end