classdef DensityAndBound < DesignVariable
    properties (Access = public)
        density
        bound
    end
    methods (Access = public)
        function obj = DensityAndBound(cParams)
            obj.density = Density(cParams);
            obj.bound = 1;
            obj.mesh = cParams.mesh;
            obj.type = cParams.type;
            obj.createValue();
        end 
        function update(obj,value)
            obj.update@DesignVariable(value);
            obj.density.value = obj.value(1:end-1);
            obj.bound = obj.value(end);
        end
        function updateOld(obj)
            obj.density.updateOld();
        end
        function v = getVariablesToPlot(obj)
            v{1} = obj.density.value;
        end       
        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.density.valFun};
            funNames = {'Density'};
        end
        function rho = computeVolumeFraction(obj)
            rho = obj.density.computeVolumeFraction(obj);
        end
        function norm = computeL2normIncrement(obj)
           norm = obj.density.computeL2normIncrement();
        end
    end
    methods (Access = private)
        function createValue(obj)
            obj.value = [obj.density.value;obj.bound];
        end
    end
end