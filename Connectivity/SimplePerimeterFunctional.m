classdef SimplePerimeterFunctional < handle

    properties (Access = private)
        mesh
        filter
        quadrature
    end

    methods (Access = public)
        function obj = SimplePerimeterFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            iter = x{2};
            x = x{1};
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            J  = obj.computeFunction(xR);
            dJ = obj.computeGradient(xR);
            dJ = obj.filter.compute(dJ,2);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.filter  = cParams.filter;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function J = computeFunction(obj, rho)
            perimeter = rho.*(1-rho);
            J      = Integrator.compute(perimeter,obj.mesh,obj.quadrature.order);
        end

        function dJ = computeGradient(obj, rho)
            dj        = 1-2*rho.fValues;
            s.fValues = dj;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            dJ        = LagrangianFunction(s);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Perimeter';
        end
    end
end