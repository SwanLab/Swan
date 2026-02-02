classdef Parameter < DesignVariable

    properties (Access = private)
        plotting
        plotter
    end

    methods (Access = public)

        function obj = Parameter(cParams)
%             obj.nVariables = 1;
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
            if size(value,2) < obj.fun.ndimf
                value = reshape(value,[],obj.fun.ndimf);
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
                % obj.plotter  = PlotterRadius(cParams);
                obj.plotter    = PlotterLatticeCross2param(cParams);
                % obj.plotter  = PlotterLatticeCross(cParams);
            end
        end

    end
    
end

