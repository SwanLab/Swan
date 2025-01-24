classdef SimplePerimeterFunctional < handle

    properties (Access = private)
        mesh
        filter
        quadrature
        weight
        updated
    end

    methods (Access = public)
        function obj = SimplePerimeterFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            iter = x{2};
            x = x{1};
            if iter == 150 && obj.updated == 0
                obj.weight = 100;
                obj.updated = 1;
            end
   
            xD = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);
            J  = obj.computeFunction(xD,xR);
            dJ = obj.computeGradient(xR);
            dJ = obj.filter.compute(dJ,2);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.weight = 0.0;
            obj.updated = 0;
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

        function J = computeFunction(obj, rho, rhoF)
            perimeter = rho.*(1-rhoF);
            J      = obj.weight*Integrator.compute(perimeter,obj.mesh,obj.quadrature.order);
        end

        function dJ = computeGradient(obj, rhoF)
            dj        = 1-2*rhoF.fValues;
            s.fValues = obj.weight*dj;
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