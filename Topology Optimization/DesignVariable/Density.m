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
            fun = obj.fun{1};
        end

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            obj.fun{1}   = obj.fun{1}.copy();
            obj.fun{1}.fValues = value;
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot(obj.fun{1});
            end
        end
    
    end

    methods (Access = private)

        function createPlotter(obj,cParams)
            obj.plotting = cParams.plotting;
            if obj.plotting
                s.type    = 'Density';
                s.mesh    = obj.fun{1}.mesh;
                obj.plotter  = Plotter.create(s);
            end
        end

    end
    
end

