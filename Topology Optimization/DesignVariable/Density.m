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
            fun = obj.fun;
        end

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            s.order   = 'P1';
            obj.fun   = LagrangianFunction(s);
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot();
            end
        end

        function rho = copy(obj)
            s.fun      = obj.fun;
            s.mesh     = obj.mesh;
            s.type     = 'Density';
            s.plotting = false;
            rho        = DesignVariable.create(s);
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

