classdef RadiusAndBound < DesignVariable

    properties (Access = public)
        radius
        bound
    end

    methods (Access = public)

        function obj = RadiusAndBound(cParams)
            obj.initBound(cParams);
            obj.createRadius(cParams);
            obj.computeValue();
        end

        function update(obj,value)
            obj.radius.update(value(1:end-1));
            obj.bound = value(end);
            obj.computeValue();
        end

        function plot(obj)
            obj.radius.plot();
        end

        function updateOld(obj)
            obj.radius.updateOld();
        end

        function norm = computeL2normIncrement(obj)
           norm = obj.radius.computeL2normIncrement();
        end

    end

    methods (Access = private)

        function initBound(obj,cParams)
            obj.bound = 1;
            obj.type  = cParams.type;
        end

        function createRadius(obj,cParams)
            cParams.type = 'Radius';
            obj.radius  = DesignVariable.create(cParams);
        end

        function computeValue(obj)
            obj.fun.fValues = [obj.radius.fun.fValues;obj.bound];
        end
        
    end
end