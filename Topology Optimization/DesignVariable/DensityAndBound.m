classdef DensityAndBound < DesignVariable
    properties (Access = public)
        value
    end
    properties (Access = private)
        density
        bound
    end
    methods (Access = public)
        function obj = DensityAndBound(cParams)
            obj.density = Density(cParams);
            obj.bound = 5;
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
    end
    methods (Access = private)
        function createValue(obj)
            obj.value = [obj.density.value;obj.bound];
        end
    end
end