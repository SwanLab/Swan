classdef DensityAndBound < DesignVariable

    properties (Access = public)
        density
        bound
    end

    methods (Access = public)

        function obj = DensityAndBound(cParams)
            obj.initBound(cParams);
            obj.createDensity(cParams);
            obj.computeValue();
        end

        function update(obj,value)
            obj.density.update(value(1:end-1));
            obj.bound = value(end);
            obj.computeValue();
        end

        function plot(obj)
            obj.density.plot();
        end

        function updateOld(obj)
            obj.density.updateOld();
        end

        function norm = computeL2normIncrement(obj)
           norm = obj.density.computeL2normIncrement();
        end

    end

    methods (Access = private)

        function initBound(obj,cParams)
            obj.bound = 1;
            obj.type  = cParams.type;
        end

        function createDensity(obj,cParams)
            cParams.type = 'Density';
            obj.density  = DesignVariable.create(cParams);
        end

        function computeValue(obj)
            obj.fun.fValues = [obj.density.fun.fValues;obj.bound];
        end
        
    end
end