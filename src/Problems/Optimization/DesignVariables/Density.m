classdef Density < DesignVariable

    properties (Access = private)
        plotting
        plotter
    end

    methods (Access = public)

        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
            obj.createPlotter(cParams);
        end

        function fun = obtainDomainFunction(obj)
            fun{1} = obj.fun;
        end

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            obj.fun.setFValues(value);
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot();
            end
        end
    
    end

    methods (Access = private)

        function createPlotter(obj,cParams)
            obj.plotting = cParams.plotting;
            if obj.plotting
                obj.plotter  = Plotter.create(obj);
            end
        end

    end
    
end

